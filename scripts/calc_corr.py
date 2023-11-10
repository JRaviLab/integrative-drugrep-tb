import pandas as pd
import numpy as np
import os
import multiprocessing
from utils.corr_funcs import *

# FIND DISEASE CONTROL SAMPLES
control_samples = []
args = pd.read_csv('data/metadata/microarray/TB_microarray_args.tsv', delimiter='\t', header=0)
for dir in ['data/metadata/microarray/', 'data/metadata/rnaseq/']: # Iterate over directories
    if dir == 'data/metadata/rnaseq/':
        args = pd.read_csv('data/metadata/rnaseq/TB_rnaseq_args.tsv', delimiter='\t', header=0)
    for i, file_name in enumerate(args['meta_path']):
        file_name = file_name.split('/')[-1]
        if 'args' in file_name:
            continue
        metadata = pd.read_csv(os.path.join(dir, file_name), delimiter='\t', header=0)
        args_info = args.loc[i]
        keyword = args_info['control_keywords']
        col_name = args_info['condition_colname']
        for j in metadata.index: # Iterate over rows in the metadata DataFrame
            location = str(metadata.loc[j, col_name])
            try:
                if (keyword in location): # Check for keywords in the sample location
                    if dir == 'data/metadata/rnaseq/':
                        control_samples.append(metadata.iloc[j, 1])
                    else:
                        control_samples.append(metadata.iloc[j, 0])
            except:
                pass
control_samples = np.array(control_samples)
control_samples = np.unique(control_samples)
del args, metadata, i, j, file_name, dir, col_name, keyword, location, args_info

# READ IN DRUG DATA
drug_data = pd.read_csv('data/expression/drug_data/level3_control_expdata_lm.tsv', delimiter='\t')
drug_data.set_index('rid', inplace=True)
landmark_genes = drug_data.index

# READ IN DISEASE DATA
dis_meta = pd.read_csv('data/annotation/Homo_sapiens.gene_info.csv', header=[0], index_col=[0])
# generate maps from Ensembl and Symbol
ma_map = dis_meta.drop_duplicates(subset='Ensembl').set_index('Ensembl')['GeneID'] # microarray
rn_map = dis_meta.drop_duplicates(subset='Symbol').set_index('Symbol')['GeneID'] # rnaseq
# inputs a data file and return an expression df with GeneIDs of landmark genes
def disease_prep(filename, ma=True):
    dis_exp = pd.read_csv(filename, delimiter='\t', low_memory=False, dtype={'GeneID': 'int32'})
    dis_exp.columns.values[0] = 'GeneID'
    dis_exp['GeneID'] = dis_exp['GeneID'].map(ma_map if ma else rn_map)
    dis_exp = dis_exp.loc[np.isin(dis_exp['GeneID'], landmark_genes)] # subset for landmark genes
    dis_exp = dis_exp.loc[:, np.isin(dis_exp.columns.values, np.append(control_samples, 'GeneID'))] # subset for control samples
    return dis_exp
dis_data = pd.DataFrame({"GeneID": landmark_genes})
for dir in ['data/expression/microarray/', 'data/expression/rnaseq/']: # Iterate over directories
    for file_name in os.listdir(dir):
        dis_data = dis_data.merge(disease_prep(dir+file_name, dir[-2]=='y'), how='left', on='GeneID')
dis_data = dis_data.set_index('GeneID').astype(float)
dis_data.columns.name = "Disease Samples"
dis_data = dis_data.T.drop_duplicates().T
del dis_meta, ma_map, rn_map, dir, file_name

# Determine the number of CPU cores to utilize
num_cores = multiprocessing.cpu_count()
print(f"Number of CPU cores: {num_cores}")

# Split drug samples into subsets
drug_sample_chunks = np.array_split(drug_data.columns, num_cores)

# Create a multiprocessing pool
pool = multiprocessing.Pool(processes=num_cores)

# Parallelize the correlation calculations
results = pool.starmap(calculate_correlations, [(drug_data[chunk], dis_data) for chunk in drug_sample_chunks])

# Combine the results into the final correlation matrix
for i, result in enumerate(results):
    if i == 0:
        final_correlation_matrix = result
    else:
        final_correlation_matrix = pd.concat([final_correlation_matrix,result])
    
final_correlation_matrix.to_hdf('final_correlation_matrix.h5', key='cor', mode='w')

# Close the pool
pool.close()
pool.join()
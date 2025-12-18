import pandas as pd
import numpy as np
import os

def read_drug_data(quantile_normalize=True, groupby_cell_line=True):
    try:
        drug_data = pd.read_csv('data/expression/drug_data/level3_control_expdata_lm.tsv', delimiter='\t')
    except:
        drug_data = pd.read_csv('../data/expression/drug_data/level3_control_expdata_lm.tsv', delimiter='\t')
    drug_data.set_index('rid', inplace=True)
    landmark_genes = drug_data.index
    drug_samples = drug_data.columns
    
    split_drug = pd.DataFrame(columns=['exp', 'cell', 'time', 'x','end'])
    for sam in drug_samples:
        split_drug.loc[len(split_drug.index)] = sam.split('_')
    split_drug = split_drug[split_drug['time'] == '24H'] # Remove 6H and 48H: 19258 -> 18182
    split_drug.sort_values(by=['exp', 'cell', 'end'], inplace=True)
    drug_data = drug_data.loc[:,split_drug['exp']+'_'+split_drug['cell']+'_'+split_drug['time']+'_'+split_drug['x']+'_'+split_drug['end']]
    split_drug.reset_index(drop=True, inplace=True)

    def row_sample(rowi, include_x=False):
        columns = ['exp', 'cell', 'time']
        if include_x:
            columns.append('x')
        columns.append('end')
        
        return '_'.join(str(split_drug.loc[rowi, col]) for col in columns)

    past = row_sample(0)
    group = [row_sample(0, include_x=True)]
    cols_to_concat = []
    col_names = [past]

    for i in range(1, len(split_drug.index)):
        cur = row_sample(i)
        if cur == past:
            group.append(row_sample(i, include_x=True))
        else:
            cols_to_concat.append(drug_data[group].mean(axis=1))
            group = [row_sample(i, include_x=True)]
            past = cur
            col_names.append(cur)
    cols_to_concat.append(drug_data[group].mean(axis=1))

    drug_data = pd.concat(cols_to_concat, axis=1)
    drug_data.columns = col_names

    if groupby_cell_line:
        # group drug samples by cell line -> drug_data_cell
        cell_lines = list(set([x.split('_')[1] for x in drug_data.columns]))
        drug_data_cell = pd.DataFrame(index=drug_data.index, columns=cell_lines)
        for cell_line in cell_lines:
            cell_line_columns = [x for x in drug_data.columns if cell_line in x]
            drug_data_cell[cell_line] = drug_data[cell_line_columns].mean(axis=1)
        drug_data = drug_data_cell

    if quantile_normalize:
        # https://github.com/ShawnLYU/Quantile_Normalize/blob/master/quantile_norm.py
        def quantileNormalize(df_input):
            df = df_input.copy()
            #compute rank
            dic = {}
            for col in df:
                dic.update({col : sorted(df[col])})
            sorted_df = pd.DataFrame(dic)
            rank = sorted_df.mean(axis = 1).tolist()
            #sort
            for col in df:
                t = np.searchsorted(np.sort(df[col]), df[col])
                df[col] = [rank[i] for i in t]
            return df
        
        drug_data = quantileNormalize(drug_data)

    return drug_data, landmark_genes

# Uses metadata to form an array of disease sample ids that are control samples
def get_disease_controls():
    control_samples = []
    for dir in ['data/metadata/microarray/', 'data/metadata/rnaseq/']:
        if dir == 'data/metadata/microarray/':
            args = pd.read_csv('data/metadata/microarray/TB_microarray_args.tsv', 
                            delimiter='\t', header=0)
        else:
            args = pd.read_csv('data/metadata/rnaseq/TB_rnaseq_args.tsv', 
                            delimiter='\t', header=0)
        for i, file_name in enumerate(args['meta_path']):
            file_name = file_name.split('/')[-1]
            if 'args' in file_name:
                continue
            metadata = pd.read_csv(os.path.join(dir, file_name), delimiter='\t', 
                                header=0)
            args_info = args.loc[i]
            keyword = args_info['control_keywords']
            col_name = args_info['condition_colname']
            for j in metadata.index:
                location = str(metadata.loc[j, col_name])
                try:
                    if (keyword in location): # Check for keywords in the sample location
                        if dir == 'data/metadata/rnaseq/':
                            control_samples.append(metadata.iloc[j, 1])
                        else:
                            control_samples.append(metadata.iloc[j, 0])
                except:
                    pass
    return np.unique(np.array(control_samples))

def read_dis_data(landmark_genes, control_samples=get_disease_controls(), rnaseq_symbol=False):
    dis_meta = pd.read_csv('data/annotation/Homo_sapiens.gene_info.csv', header=[0], index_col=[0])
    dis_meta = dis_meta[[element in landmark_genes for element in dis_meta['GeneID'].values]]
    # generate maps from Ensembl and Symbol
    ma_map = dis_meta.drop_duplicates(subset='Ensembl').set_index('Ensembl')['GeneID'] # microarray
    if rnaseq_symbol:
        rn_map = dis_meta.drop_duplicates(subset='Symbol').set_index('Symbol')['GeneID'] # rnaseq
    else:
        rn_map = ma_map
    # inputs a data file and return an expression df with GeneIDs of landmark genes
    def disease_prep(filename, ma=True):
        dis_exp = pd.read_csv(filename, delimiter='\t', low_memory=False, dtype={'GeneID': 'int32'})
        dis_exp.columns.values[0] = 'GeneID'
        dis_exp['GeneID'] = dis_exp['GeneID'].map(ma_map if ma else rn_map)
        dis_exp = dis_exp.loc[np.isin(dis_exp['GeneID'], landmark_genes)] # subset for landmark genes
        dis_exp = dis_exp.loc[:, np.isin(dis_exp.columns.values, np.append(control_samples, 'GeneID'))] # subset for control samples
        return dis_exp
    rnaseq_dis_data = pd.DataFrame({"GeneID": landmark_genes})
  
    for file_name in os.listdir('data/expression/rnaseq/'):
        rnaseq_dis_data = rnaseq_dis_data.merge(disease_prep('data/expression/rnaseq/'+file_name, ma=False), how='left', on='GeneID')
    rnaseq_dis_data = rnaseq_dis_data.set_index('GeneID').astype(float)
    rnaseq_dis_data.columns.name = "Disease Samples"
    rnaseq_dis_data = rnaseq_dis_data.T.drop_duplicates().T
    
    marray_dis_data = pd.DataFrame({"GeneID": landmark_genes})
    for file_name in os.listdir('data/expression/microarray/'):
        marray_dis_data = marray_dis_data.merge(disease_prep('data/expression/microarray/'+file_name, ma=True), how='left', on='GeneID')
    marray_dis_data = marray_dis_data.set_index('GeneID').astype(float)
    marray_dis_data.columns.name = "Disease Samples"
    marray_dis_data = marray_dis_data.T.drop_duplicates().T
    return rnaseq_dis_data, marray_dis_data

def get_dis_meta(dis_sample, label):
    dis_meta = pd.read_csv('dis_meta.csv', sep=',', header=0)
    return dis_meta.loc[dis_meta['dis_sample'] == dis_sample, label].values[0]
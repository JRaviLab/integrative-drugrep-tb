# script to extract expression data of all untreated samples in LINCS data
# last modified: 03/13/26
# Kewalin Samart

from cmapPy.pandasGEXpress import parse
import pandas as pd
import numpy as np
import shelve

# read in metadata 
lincs_level3_meta = pd.read_csv("../data/metadata/GSE92742_Broad_LINCS_inst_info.txt",sep="\t")
# read in landmark genes
landmark_genes = list(pd.read_csv("../data/metadata/LINCS_landmark_gene_annotation.csv")["pr_gene_id"])
landmark_genes = [str(gene_id) for gene_id in landmark_genes]

# grab all inst_id with pert_type == "ctl_vehicle"
control_df = lincs_level3_meta[lincs_level3_meta["pert_type"] == "ctl_vehicle"] 
control_samples = list(control_df["inst_id"])

# parsing gctx file to gct object to get expression matrix
level3_gct = parse.parse("../data/baseline_analysis/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx")
# level3_gct contains 3 data components: (i) row_metadata_df (genes), (ii) column_metadata_df (samples), and (iii) data_df (expression data)
level3_control_expdata = level3_gct.data_df[control_samples]
# subset data to include only landmark genes
level3_control_expdata_lm = level3_control_expdata.loc[landmark_genes]

# write out the dataframe of control samples for drug data
level3_control_expdata_lm.to_csv('../data/baseline_analysis/level3_control_expdata_lm.tsv', sep='\t', index=True)

print(level3_control_expdata_lm.columns)
print(level3_control_expdata_lm.index)

# saving row (genes; entrezid) and column (sample ids) names 
with shelve.open('../data/baseline_analysis/level3_drugdata_rows_cols_info') as shelveFile:
    shelveFile['sample_ids'] = level3_control_expdata_lm.columns # 19258 samples
    shelveFile['genes'] = level3_control_expdata_lm.index # 978 landmark genes
    shelveFile.close()

# numpy array --> condense to small file 
np_level3_control_expdata_lm = level3_control_expdata_lm.to_numpy()

np.save('../data/baseline_analysis/level3_control_expdata_lm.npy', np_level3_control_expdata_lm, allow_pickle=True, fix_imports=True)
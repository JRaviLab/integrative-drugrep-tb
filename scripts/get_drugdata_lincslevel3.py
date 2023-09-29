# script to extract expression data of all untreated samples in LINCS data
# created date: 09/29/23
# Kewalin Samart

from cmapPy.pandasGEXpress import parse
import pandas as pd

# read in metadata 
lincs_level3_meta = pd.read_csv("../data/database/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt",sep="\t")
# read in landmark genes
landmark_genes = list(pd.read_csv("../data/annotation/LINCS_landmark_gene_annotation.csv")["pr_gene_id"])
landmark_genes = [str(gene_id) for gene_id in landmark_genes]

# grab all inst_id with pert_type == "ctl_vehicle"
control_df = lincs_level3_meta[lincs_level3_meta["pert_type"] == "ctl_vehicle"] 
control_samples = list(control_df["inst_id"])

# parsing gctx file to gct object to get expression matrix
level3_gct = parse.parse("../data/database/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx")
# level3_gct contains 3 data components: (i) row_metadata_df (genes), (ii) column_metadata_df (samples), and (iii) data_df (expression data)
level3_control_expdata = level3_gct.data_df[control_samples]
# subset data to include only landmark genes
level3_control_expdata_lm = level3_control_expdata.loc[landmark_genes]

# write out the dataframe of control samples for drug data
level3_control_expdata_lm.to_csv('../data/expression/drug_data/level3_control_expdata_lm.tsv', sep='\t', index=True)


# script to generate data matrix in h5 format via CmapPy
# created date: 01/23/24
# last modified: 03/08/24
# Kewalin Samart

library(signatureSearch)
library(cmapR)
library(readr)


# read in gctx LINCS data level 5
# path to gctx file
gctx_path = "/projects/ksamart@xsede.org/drugrep_tb/data/drug_data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gct <- parse_gctx(gctx_path)
#print(head(gct))
#save.image('session_w_gct.rda')
subset_gct_test <- cmapR::subset_gct(gct,cid=1)
drug_level5_mat <- gct@mat
drug_level5_rdesc <- gct@rdesc
drug_level5_cdesc <- gct@cdesc

#save.image('session_LINCSdruglevel5.rda')

#write_tsv(as.data.frame(drug_level5_mat), file = "/data/drug_data/drug_level5_mat.tsv")

# read in L1000 GeneIDs
L1000_data <- read_csv("/projects/ksamart@xsede.org/drugrep_tb/data/metadata/LINCS_gene_annotation.csv")
L1000_genes <- L1000_data$pr_gene_id

# subset only L1000 genes
level5_data_landmark = drug_level5_mat[rownames(drug_level5_mat) %in% L1000_genes,] # filtering out non L1000 genes

# convert data to .h5 format
# reformat column names to (drug)__(cell line)__(conditions)
colnames(level5_data_landmark) <- gsub("_","__",colnames(level5_data_landmark))
h5file <- "/projects/ksamart@xsede.org/drugrep_tb/data/drug_data/level5_data_landmark.h5"
# save .h5 file
build_custom_db(level5_data_landmark, h5file)

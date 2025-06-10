# script to finalize enriched parent terms from the summarized pathway results
# add disease-drug pathway-associated genes to the drug-target network
# created date: 04/15/24
# created date: 04/16/24
# Kewalin Samart

library(readr)
library(rrvgo)
library(GO.db)
library(biomaRt)
library(org.Hs.eg.db)

# read in representative GO matrices -- microarray and RNAseq of the same regulation direction (up/dn)
rep_GO_marray_up <- as.data.frame(read_tsv("./data/pathways/microarray/up/GO_ORA/rep_GO_mat_up.tsv"))
rep_GO_rnaseq_up <- as.data.frame(read_tsv("./data/pathways/RNAseq/up/GO_ORA/rep_GO_mat_up.tsv"))

rep_GO_marray_dn <- as.data.frame(read_tsv("./data/pathways/microarray/dn/GO_ORA/rep_GO_mat_dn.tsv"))
rep_GO_rnaseq_dn <- as.data.frame(read_tsv("./data/pathways/RNAseq/dn/GO_ORA/rep_GO_mat_dn.tsv"))

# only have numerical values in the entries
# upregulated pathways
row.names(rep_GO_marray_up) <- rep_GO_marray_up$ID
rep_GO_marray_up$Description <- NULL
rep_GO_marray_up$ID <- NULL
# get pathways apperaing in at least 40% of the signatures
rep_enrichedGO_marray_up <- rep_GO_marray_up[rowSums(rep_GO_marray_up == 0) <= 0.8*ncol(rep_GO_marray_up),]

row.names(rep_GO_rnaseq_up) <- rep_GO_rnaseq_up$ID
rep_GO_rnaseq_up$Description <- NULL
rep_GO_rnaseq_up$ID <- NULL
# get pathways apperaing in at least 20% of the signatures
rep_enrichedGO_rnaseq_up <- rep_GO_rnaseq_up[rowSums(rep_GO_rnaseq_up == 0) <= 0.8*ncol(rep_GO_rnaseq_up),]

# downregulated pathways
row.names(rep_GO_marray_dn) <- rep_GO_marray_dn$ID
rep_GO_marray_dn$Description <- NULL
rep_GO_marray_dn$ID <- NULL
# get pathways apperaing in at least 20% of the signatures
rep_enrichedGO_marray_dn <- rep_GO_marray_dn[rowSums(rep_GO_marray_dn == 0) <= 0.8*ncol(rep_GO_marray_dn),]

row.names(rep_GO_rnaseq_dn) <- rep_GO_rnaseq_dn$ID
rep_GO_rnaseq_dn$Description <- NULL
rep_GO_rnaseq_dn$ID <- NULL
# get pathways apperaing in at least 20% of the signatures
rep_enrichedGO_rnaseq_dn <- rep_GO_rnaseq_dn[rowSums(rep_GO_rnaseq_dn == 0) <= 0.8*ncol(rep_GO_rnaseq_dn),]

# try intersection of both microarray and RNAseq and see how many terms we have left
# if not much, then use the union instead
rep_GO_up <- union(row.names(rep_GO_marray_up),row.names(rep_GO_rnaseq_up))
rep_GO_dn <- union(row.names(rep_GO_marray_dn),row.names(rep_GO_rnaseq_dn))

rep_enrichedGO_up <- union(row.names(rep_enrichedGO_marray_up),row.names(rep_enrichedGO_rnaseq_up))
rep_enrichedGO_dn <- union(row.names(rep_enrichedGO_marray_dn),row.names(rep_enrichedGO_rnaseq_dn))

# go back to the genes contributing to these pathways starting with the genes with the highest number of occurrence
rep_GOterm_up <- getGOTerm(rep_GO_up)$BP
rep_GOterm_dn <- getGOTerm(rep_GO_dn)$BP

rep_enrichedGOterm_up <- getGOTerm(rep_enrichedGO_up)
rep_enrichedGOterm_dn <- getGOTerm(rep_enrichedGO_dn)

# read in DrugTarget_nodes_unweighted.tsv
DrugTarget_network <- read_tsv("./results/uniformly_processed/DrugTarget_nodes_unweighted.tsv")
drug_pathways <- DrugTarget_network$pathway

drug_dis_pathways_up <- intersect(rep_GOterm_up, drug_pathways) # GO:0019373
drug_dis_pathways_dn <- intersect(rep_GOterm_dn, drug_pathways) # GO:0019373, GO:0070988

# get children terms
getGOChildren("GO:0019373")
getGOChildren(c("GO:0019373", "GO:0070988"))$`GO:0070988`$Children

mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
symb <- keys(org.Hs.eg.db, "SYMBOL")
genes.0019373_goid <- getBM(attributes=c('hgnc_symbol', 'go_id'),values = 'GO:0019373', mart = mart)
genes.0070988_goid <- getBM(attributes=c('hgnc_symbol', 'go_id'),values = 'GO:0070988', mart = mart)

gene.0019373_vec <- genes.0019373_goid[genes.0019373_goid$go_id == 'GO:0019373',]$hgnc_symbol
gene.0070988_vec <- genes.0070988_goid[genes.0070988_goid$go_id == 'GO:0070988',]$hgnc_symbol

drug_dis_genes <- c(gene.0019373_vec,gene.0070988_vec)

# overlay the genes on the STRING network
# read in STRING interaction network
STRING_path = "./data/metadata/STRING_hs_interactions.txt"
gene_gene_string <- read.delim(STRING_path,col.names = c("gene1","gene2","edge_weight"))
# get protein-protein interactions
gene_gene_string_subset1 <- gene_gene_string[which(gene_gene_string$gene1 %in% drug_dis_genes),]
gene_gene_string_finalsubset <- gene_gene_string_subset1[which(gene_gene_string_subset1$gene2 %in% drug_dis_genes),]
# get unique interactions
gene_gene_string_finalsubset <- gene_gene_string_finalsubset[!duplicated(gene_gene_string_finalsubset$edge_weight),]
gene_gene_string_finalsubset$pathway = "epoxygenase P450 pathway"
gene_gene_string_finalsubset <- gene_gene_string_finalsubset[c("pathway","gene1","gene2")]
colnames(gene_gene_string_finalsubset)[2] <- "drugdis_gene1"
colnames(gene_gene_string_finalsubset)[3] <- "drugdis_gene2"
# merge with DrugTarget_network using pathways
colnames(DrugTarget_network) <- gsub("g","moa_g",colnames(DrugTarget_network))
DrugTarget_network_dis <- merge(DrugTarget_network,gene_gene_string_finalsubset,by="pathway",all.x=TRUE,all.y=TRUE)
write_tsv(DrugTarget_network_dis,file="./results/uniformly_processed/DrugTarget_network_disgenes-pathway.tsv")

drugdis_edgeinfo <- gene_gene_string_finalsubset[c("drugdis_gene1","drugdis_gene2")]
drugdis_edgeinfo$type <- "dd_gene-dd_gene"
colnames(drugdis_edgeinfo) <- c("from","to","type")

pathway_gene1_edgeinfo <- gene_gene_string_finalsubset[c("pathway","drugdis_gene1")]
pathway_gene1_edgeinfo$type <- "pathway-dd_gene"
colnames(pathway_gene1_edgeinfo) <- c("from","to","type")

pathway_gene2_edgeinfo <- gene_gene_string_finalsubset[c("pathway","drugdis_gene2")]
pathway_gene2_edgeinfo$type <- "pathway-dd_gene"
colnames(pathway_gene2_edgeinfo) <- c("from","to","type")

DrugDisTarget_network_edgeinfo <- rbind(drugdis_edgeinfo,pathway_gene1_edgeinfo)
DrugDisTarget_network_edgeinfo <- rbind(DrugDisTarget_network_edgeinfo,pathway_gene2_edgeinfo)
DrugDisTarget_network_edgeinfo <- unique(DrugDisTarget_network_edgeinfo)

# read in drug-target edge info
drugtarget_edges <- read_tsv("./results/uniformly_processed/edges_unweighted.tsv")
DrugDisTarget_network_edgeinfo <- rbind(drugtarget_edges, DrugDisTarget_network_edgeinfo)

write_tsv(DrugDisTarget_network_edgeinfo,file="./results/uniformly_processed/DrugDisTarget_network_edgeinfo.tsv")

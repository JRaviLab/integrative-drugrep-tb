# script to perform pathway analyses: GO/KEGG/Reactome ORA
# created date: 06/25/22
# last modified: 04/24/24
# Kewalin Samart

library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(readr)

source("./scripts/00_background_genes_PA_functions.R")

# set up arguments
args <- commandArgs(TRUE)
path_to_sigs <- args[1]
direction <- args[2]
pattern <- args[3]
technology <- args[4]
bg_genes_source <- args[5]
extra_arg <- args[6]
metadata_path <- args[7]

## example arguments
#path_to_sigs <- "./data/uniformly_processed/RNAseq/signatures"
#direction <- "dn" #"dn", "full" or ""
#pattern <- "dn_none.tsv"
#technology <- "RNAseq" #"microarray"
#bg_genes_source <- "input data" # "LINCS", "KEGG", "GO", "input data"
#extra_arg <- "" # a string or a list of strings indicating one or more types of LINCS gene: (i) landmark (by default) (ii) inferred (iii) best inferred (iv) not inferred (v) reference
#metadata_path <- "./inputs/TB_rnaseq_args.tsv"
#extension <- "dn_none"


# Microarray up -- done # rerun for old bg genes from input data only
# Microarray dn -- done
# RNAseq up -- done
# RNAseq dn

# path to DE data table/ signatures
data_path <- paste0(path_to_sigs,"/",direction)

# get file names
filenames <- list.files(path = data_path,
                        pattern = pattern, all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

dirname__ <- paste0("./data/pathways/")
if(!dir.exists(dirname__)){
  dir.create(dirname__)
}

dirname_ <- paste0("./data/pathways/",technology)
if(!dir.exists(dirname_)){
  dir.create(dirname_)
}

dirname <- paste0("./data/pathways/",technology,"/",direction)
if(!dir.exists(dirname)){
  dir.create(dirname)
}

dirname_read <- paste0(path_to_sigs,"/",direction)

if(direction == "up"){
  scoreType <- "pos"
}else if(direction == "dn"){
  scoreType <- "neg"
}else if(direction == "full"){
  scoreType <- "std"
}

# define background genes to use
if(extra_arg == ""){extra_arg = NULL}
bg_genes_ <- get_bg_genes(bg_source=bg_genes_source, metadata_path=metadata_path, data_path=data_path, extra_arg=extra_arg, extension=extension)
bg_genes <- intersect(bg_genes_,GO_genes()) # intersect with GO genes

itr = 1
for(file in filenames){

  print(paste0("Running iteration: ",itr))

  genes_df <- read_tsv(paste0(dirname_read,"/",file))
  if(direction == "full"){
    gene_vector <- genes_df$log2FoldChange
    names(gene_vector) <- as.character(genes_df$GeneID)
  }else if(direction == ""){
    gene_vector <- genes_df$aggregated_GeneScores
    names(gene_vector) <- as.character(genes_df$GeneID)
  }else{
    gene_vector <- genes_df$log2FoldChange
    names(gene_vector) <- as.character(genes_df$GeneID)
  }
  ##############
  # ORA GO
  print(paste0("Getting enriched GO pathways for ", file))

  enrichGO_res <- enrichGO(gene = names(gene_vector),
                           OrgDb    = 'org.Hs.eg.db',
                           readable = T,
                           ont = "BP",
                           universe = bg_genes,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.1,
                           minGSSize = 5,
                           maxGSSize = 200)
  enrichGO_res <- as.data.frame(enrichGO_res)

  dirname_GO <- paste0(dirname,"/GO_ORA")
  if(!dir.exists(dirname_GO)){
    dir.create(dirname_GO)
  }

  write_tsv(enrichGO_res, file = paste0(dirname_GO,"/GO_ORA_BGcorrected_",file))
  ############## ** comment out KEGG and Reactome as we want just GO terms
  # ORA KEGG
  #print(paste0("Getting enriched KEGG pathways for ", file))

  #enrichKEGG_res <- enrichKEGG(gene = names(gene_vector),
  #                             organism = 'hsa',
  #                             universe = bg_genes,
  #                             pvalueCutoff = 0.05,
  #                             qvalueCutoff = 0.1,
  #                             minGSSize = 5,
  #                             maxGSSize = 200)
  #enrichKEGG_res <- as.data.frame(enrichKEGG_res)

  #dirname_KEGG <- paste0(dirname,"/KEGG_ORA")
  #if(!dir.exists(dirname_KEGG)){
  #  dir.create(dirname_KEGG)
  #}

  #write_tsv(enrichKEGG_res, file = paste0(dirname_KEGG,"/KEGG_ORA_",file))
  #############
  # ORA Reactome
  #print(paste0("Getting enriched Reactome pathways for ", file))

  #enrichReactome_res <- enrichPathway(gene = names(gene_vector),
  #                                    organism = "human",
  #                                    pvalueCutoff = 0.05,
  #                                    universe = bg_genes,
  #                                    readable=T,
  #                                    minGSSize = 5,
  #                                    maxGSSize = 200)
  #enrichReactome_res <- as.data.frame(enrichReactome_res)

  #dirname_Reactome <- paste0(dirname,"/Reactome_ORA")
  #if(!dir.exists(dirname_Reactome)){
  #  dir.create(dirname_Reactome)
  #}

  #write_tsv(enrichReactome_res, file = paste0(dirname_Reactome,"/Reactome_ORA_",file))

  itr = itr + 1
}


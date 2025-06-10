# script to get a shortlist of pathway representatives for a given GO results
# created date: 07/25/22
# last modified: 04/15/24
# Kewalin Samart

library(rrvgo)
library(pheatmap)
library(readr)

# set up arguments
args <- commandArgs(TRUE)
data_path <- args[1] # e.g. "./data/pathways/RNAseq/dn/GO_ORA"
direction <- args[2] # "dn"
pattern <- args[3] # "dn_none" "up_none_correctedBG.tsv"

# get file names
filenames <- list.files(path = data_path,
                        pattern = pattern, all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
union_terms <- c()
for(i in 1:length(filenames)){

  file = filenames[i]
  GO_res <- read_tsv(paste0(data_path,"/",file))

  # get pathway representatives
  print("Getting pathway representatives")
  print(paste0("iteration ",i,": ", file))

  tryCatch({
    simMatrix <- calculateSimMatrix(GO_res$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
    scores <- setNames(-log10(GO_res$qvalue), GO_res$ID)
    reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
  }, error=function(e){})

  write_tsv(reducedTerms,file = paste0(data_path,"/summarized_pathways/reducedTerms_",file))

  reduced_GO <- GO_res[GO_res$ID %in% unique(reducedTerms$parent),]

  write_tsv(reduced_GO,file = paste0(data_path,"/summarized_pathways/parentTerms_",file))

  # get only data description i.e., signature name
  signature_name <- gsub("GO_ORA_","",file)
  signature_name <- gsub(paste0("_",pattern,".tsv"),"",signature_name)

  reduced_GO_col <- reduced_GO[c("ID","qvalue")] # debug starts here # look for a mapping file GO terms <-> description
  colnames(reduced_GO_col)[2] <- signature_name

  if(i == 1){
    rep_GO_mat <- reduced_GO_col
    union_terms <- reduced_GO$ID
  }else{
    rep_GO_mat <- merge(rep_GO_mat,reduced_GO_col, by = "ID", all.x = TRUE, all.y = TRUE)
    union_terms <- union(union_terms,reduced_GO$ID)
  }
  print(paste0("Numer of unique pathway representative across ",i," signatures: ",length(union_terms)))
}

# finalize matrix
rep_GO_mat_ <- rep_GO_mat
rep_GO_mat_$ID <- NULL
rep_GO_mat_[is.na(rep_GO_mat_)] <- 0.0
rep_GO_mat_ <- as.matrix(rep_GO_mat_)
rep_GO_mat_ <- -log10(rep_GO_mat_)
rep_GO_mat_[is.infinite(rep_GO_mat_)] <- 0.0

# add GO description
# use getGoTerm to annotate GO terms to description
rep_GO_mat_df <- as.data.frame(rep_GO_mat_)
rep_GO_mat_df$ID <- rep_GO_mat$ID
GO_annotation <- getGOTerm(rep_GO_mat_df$ID)$BP
rep_GO_mat_df$Description <- GO_annotation
rep_GO_mat_df <- rep_GO_mat_df[c(tail(colnames(rep_GO_mat_df),2),head(colnames(rep_GO_mat_df),length(filenames)))]

write_tsv(rep_GO_mat_df, file = paste0(data_path,"/rep_GO_mat_",direction,".tsv"))

# script to get drug result summary and top drugs from an aggregated signature
# created date: 07/20/22
# last modified: 09/02/24
# Kewalin Samart

library(readr)

source("./scripts/03_summarize_drugs_methodswise_functions.R")
# set up arguments
args <- commandArgs(TRUE)

drug_summ_args <- args[1] # "./inputs/aggrSig_drug_summarization_rnaseq_args.tsv"
prefix <- args[2] # "uni_rnaseq"
output_dir <- args[3] # "./results/uniformly_processed/RNAseq/"

# read in metadata for drug results derived from aggregated disease signature
data_to_run <- read.delim(drug_summ_args, sep="\t")
print("------------------ Drug Summarization of Aggregated-Signatures Results ------------------")
print(paste0("Number of drug results to summarize: ",dim(data_to_run)[1]))

# set up output directory
# make the directory to write the files to (needs to be high I/O capable)
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

all_signi_drugs <- c()
signi_drugs_score_list <- list()

for(i in 1:nrow(data_to_run)){
    drug_res_path <- data_to_run$drug_res_path[i]
    score_method <- data_to_run$score_method[i]
    score <- data_to_run$score[i]
    print(score)
    values <- data_to_run$values[i]
    tail <- data_to_run$tail[i]
    score_percentile <- data_to_run$score_percentile[i]

    # read in drug result table
    drug_df <- read.delim(drug_res_path,sep="\t")

    # grab score vector from the result tables
    if(score == "Tau"){
        score_vec <- drug_df$Tau
        threshold <- determine_threshold(score_vec=score_vec, values=values, tail=tail, score_percentile=score_percentile)
        if(values == 'neg'){
          signi_drug_df <- drug_df[which(drug_df$Tau < threshold),]
        }else if(values == 'pos'){
          signi_drug_df <- drug_df[which(drug_df$Tau > threshold),]
        }else if(values == 'all'){
          signi_drug_df <- drug_df[which(abs(drug_df$Tau) < abs(threshold)),]
        }
        signi_drug_df <- signi_drug_df[c('pert','cell','Tau')]
    }
    if(score == "NCS"){
        score_vec <- drug_df$NCS
        threshold <- determine_threshold(score_vec=score_vec, values=values, tail=tail, score_percentile=score_percentile)
        if(values == 'neg'){
          signi_drug_df <- drug_df[which(drug_df$NCS < threshold),]
        }else if(values == 'pos'){
          signi_drug_df <- drug_df[which(drug_df$NCS > threshold),]
        }else if(values == 'all'){
          signi_drug_df <- drug_df[which(abs(drug_df$NCS) < abs(threshold)),]
        }
        signi_drug_df <- signi_drug_df[c('pert','cell','NCS')]
    }
    if(score == "WCS"){
        score_vec <- drug_df$WTCS
        threshold <- determine_threshold(score_vec=score_vec, values=values, tail=tail, score_percentile=score_percentile)
        if(values == 'neg'){
          signi_drug_df <- drug_df[which(drug_df$WTCS < threshold),]
        }else if(values == 'pos'){
          signi_drug_df <- drug_df[which(drug_df$WTCS > threshold),]
        }else if(values == 'all'){
          signi_drug_df <- drug_df[which(abs(drug_df$WTCS) < abs(threshold)),]
        }
        signi_drug_df <- signi_drug_df[c('pert','cell','WTCS')]
    }
    if(score == "CMAP"){
        score_vec <- drug_df$scaled_score
        threshold <- determine_threshold(score_vec=score_vec, values=values, tail=tail, score_percentile=score_percentile)
        if(values == 'neg'){
          signi_drug_df <- drug_df[which(drug_df$scaled_score < threshold),]
        }else if(values == 'pos'){
          signi_drug_df <- drug_df[which(drug_df$scaled_score > threshold),]
        }else if(values == 'all'){
          signi_drug_df <- drug_df[which(abs(drug_df$scaled_score) < abs(threshold)),]
        }
        signi_drug_df <- signi_drug_df[c('pert','cell','scaled_score')]
    }
    if(score == "Cor_spearman" | score == "Cor_pearson"){
        score_vec <- drug_df$cor_score
        threshold <- determine_threshold(score_vec=score_vec, values=values, tail=tail, score_percentile=score_percentile)
        if(values == 'neg'){
          signi_drug_df <- drug_df[which(drug_df$cor_score < threshold),]
        }else if(values == 'pos'){
          signi_drug_df <- drug_df[which(drug_df$cor_score > threshold),]
        }else if(values == 'all'){
          signi_drug_df <- drug_df[which(abs(drug_df$cor_score) < abs(threshold)),]
        }
        signi_drug_df <- signi_drug_df[c('pert','cell','cor_score')]

        if(score == "Cor_spearman"){
          colnames(signi_drug_df)[3] <- "spearman"
        }else if(score == "Cor_pearson"){
          colnames(signi_drug_df)[3] <- "pearson"
        }
    }
    print(paste0("threshold: ", threshold))
    if(i == 1){
      signi_drugs_score_df <- signi_drug_df
    }else{
      signi_drugs_score_df <- merge(signi_drugs_score_df, signi_drug_df, by = c('pert','cell'), all.x = TRUE, all.y = TRUE)
    }
    signi_drugs_score_list[[score]] <- signi_drug_df$pert # named list containing significant drugs prioritized by different metrics
    all_signi_drugs <- c(all_signi_drugs,signi_drug_df$pert)
}

# write drug-cell and their connectivity scores
write_tsv(signi_drugs_score_df,file = paste0(output_dir,prefix,"_aggrSig_drugs_scores_",values,"_",tail,"_",score_percentile,".tsv"))

# get all significant drugs
all_signi_drugs_df <- data.frame(unique(all_signi_drugs))
colnames(all_signi_drugs_df)[1] <- "significant_drug"

print("Summarizing drug results based on their occurrence of being prioritized by 3 methods: CMAP, LINCS, and Cor methods")

for(score in names(signi_drugs_score_list)){
  score_signi_drugs_boolean <- all_signi_drugs_df$significant_drug %in% signi_drugs_score_list[[score]]
  score_signi_drugs_numeric <- replace(score_signi_drugs_boolean, score_signi_drugs_boolean==TRUE, 1)
  all_signi_drugs_df[score] <- score_signi_drugs_numeric
}

# summarize based on 3 methods: CMAP, LINCS, and Cor methods
all_signi_drugs_df$LINCS = rowSums(all_signi_drugs_df[,c("WCS", "NCS", "Tau")])
all_signi_drugs_df$Cor = rowSums(all_signi_drugs_df[,c("Cor_spearman", "Cor_pearson")])
all_signi_drugs_df$occurrence = rowSums(all_signi_drugs_df[,c("CMAP", "LINCS", "Cor")])
signi_drug_summary_df <- all_signi_drugs_df[,c("significant_drug","CMAP","LINCS","Cor","occurrence")]

drug_info_df <- read.csv(file="./data/metadata/repurposing_drugs_20200324.csv",skip=9)
merged_signi_drug_summary <- merge(signi_drug_summary_df,
                                   drug_info_df,
                                   by.x = "significant_drug",
                                   by.y = "pert_iname",
                                   all.x = TRUE,
                                   all.y = FALSE)
merged_signi_drug_summary <- merged_signi_drug_summary[order(-merged_signi_drug_summary$occurrence),]
rownames(merged_signi_drug_summary) <- 1:nrow(merged_signi_drug_summary)

write_tsv(merged_signi_drug_summary ,file = paste0(output_dir,prefix,"_aggr_drug_summary.tsv"))

# replace non-zero values with 1
aggr_drugs_score <- merged_signi_drug_summary[,c("CMAP","LINCS","Cor")]
aggr_drugs_score[-1] <- as.integer(aggr_drugs_score[-1] != 0)
# subset only drugs that show up in at least 2 out of the 3 methods
aggr_drugs_score$medthod_freq <- rowSums(aggr_drugs_score)
aggr_drugs_res_top <- merged_signi_drug_summary[which(aggr_drugs_score$medthod_freq >= 2),]
rownames(aggr_drugs_res_top) <- 1:nrow(aggr_drugs_res_top)

write_tsv(aggr_drugs_res_top,file = paste0(output_dir,prefix,"_aggr_top_drugs.tsv"))


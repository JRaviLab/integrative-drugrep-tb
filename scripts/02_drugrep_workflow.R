# script to compute connectivity scores and get drug results given input individual disease signatures
# method choices: "LINCS", "CMAP", "Cor_spearman", "Cor_pearson"
# created date: 10/05/23
# last modified: 11/27/23
# Kewalin Samart

# import needed functions
source("./scripts/01_signature_aggregation_functions.R")
source("./scripts/01_signatureSearch_connectivity_scores_functions.R")

# set up arguments
args <- commandArgs(TRUE)
sig_metadata_path <- args[1] # e.g. "./inputs/TB_microarray_args.tsv"
sig_data_path <- args[2] # e.g., "./data/uniformly_processed/microarray/signatures/"
drugdb_name <- args[3] # "LINCS", "CMAP"
score_method <- args[4] # "LINCS", "CMAP", "Cor_spearman", "Cor_pearson"
output_dir <- args[5] # e.g., "./results/uniformly_processed/microarray/LINCS/"
bg_source <- args[6] # name of the source for background genes to use: "LINCS", "KEGG", "GO", "input data"
platform <- args[7] # "microarray" or "RNAseq"
threshold <- args[8] # aggregated gene score cutoff; 0.4 by default
extra_arg <- args[9] # extra argument; optional if "LINCS" is specified for bg_source; could be one of the options below or a combination of them as a single string with comma:
                     # (i) "landmark" (by default) (ii) "inferred" (iii) "best inferred" (iv) "not inferred" (v) "reference" for examples: "landmark" or "landmark,inffered,best inferred"
                     # if "input data" is specified for bg_source, this arg could be one of the followings: "up", "dn", "full", and "" (by default)

# read in metadata file
data_to_run <- read.delim(sig_metadata_path, sep="\t")
data_to_run <- data_to_run[which(data_to_run$signature == 1), ]
print(paste0("Number of signtures to aggregate: ",dim(data_to_run)[1]))

# set up drug database
if(drugdb_name == "CMAP"){
    db_path <- setup_cmap1db()
}else if(drugdb_name %in% c("LINCS")) {
    db_path <- setup_lincsdb()
}

## --------- aggregated signature ----------- ##
# set threshold to 0.4 if not specified
if(is.null(threshold)){
    threshold <- 0.4
}

# compute aggregate signatures for all the direction: "up", "dn", and
# concatenate into "full" aggregated signature
directions = c("up","dn")
aggr_sig_list = list()
i = 1
for(direction in directions){
    print(paste0("Computing aggregated ", direction, "signature from ", sig_metadata_path))
    # modify sig_data_path to include direction
    sig_data_direction_path <- paste0("./data/uniformly_processed/",platform,"/signatures/",direction,"/")
    membership_matrix = compute_membership_matrix(data_to_run = data_to_run, data_path=sig_data_direction_path, direction, bg_source, output_dir=sig_data_direction_path, extra_arg)
    jaccard_matrix = compute_jaccard_matrix(data_to_run = data_to_run, data_path=sig_data_direction_path,direction, output_dir=sig_data_direction_path)
    aggr_signature = aggregate_signatures(membership_matrix, jaccard_matrix, output_dir=sig_data_direction_path, threshold=threshold)
    aggr_sig_list[[i]] <- aggr_signature
    i = i + 1
}

# extract aggr signatures from the list
up_aggr_sig_df <- as.data.frame(aggr_sig_list[[1]])
dn_aggr_sig_df <- as.data.frame(aggr_sig_list[[2]])

# sort genes by their aggregated gene scores
up_aggr_sig_df <- up_aggr_sig_df[order(up_aggr_sig_df$aggregated_GeneScores, decreasing = TRUE),]
dn_aggr_sig_df <- dn_aggr_sig_df[order(dn_aggr_sig_df$aggregated_GeneScores, decreasing = FALSE),]

up_aggr_genes <- up_aggr_sig_df$GeneID
dn_aggr_genes <- dn_aggr_sig_df$GeneID

# set up output directory
# make the directory to write the files to (needs to be high I/O capable)
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# quantify drug candidates for the aggregated signature
if(score_method == "CMAP"){
    final_aggr_res <- compute_CMap1_scores(up_aggr_genes, dn_aggr_genes, db_path)
}else if (score_method == "LINCS"){
    final_aggr_res <- compute_CMap2lincs_scores(up_aggr_genes, dn_aggr_genes, db_path)
}else if(score_method %in% c("Cor_spearman","Cor_pearson")){
    # concatenate up amd dn aggregated signatures into a full signature
    full_aggr_signature <- rbind(up_aggr_sig_df, dn_aggr_sig_df)
    # deal with duplicates
    dup_genes <- full_aggr_signature$GeneID[duplicated(full_aggr_signature$GeneID)]
    for(gene in dup_genes){
        aggr_gene_scores <- full_aggr_signature$aggregated_GeneScores[full_aggr_signature$GeneID == gene]
        new_aggr_gene_score <- mean(aggr_gene_scores)
        full_aggr_signature <- full_aggr_signature[!(full_aggr_signature$GeneID == gene),]
        new_val_df = data.frame(aggregated_GeneScores = new_aggr_gene_score, GeneID=gene)
        row.names(new_val_df) <- gene
       full_aggr_signature = rbind(full_aggr_signature,new_val_df)
    }
    # sort genes in descending order by the aggregated_GeneScores
    full_aggr_signature <- full_aggr_signature[order(full_aggr_signature$aggregated_GeneScores, decreasing = TRUE), ]
    # save the computed aggregated full signature
    saveRDS(full_aggr_signature, file = paste0(sig_data_path,"full/","full_aggregated_signature.rds"))
    write_tsv(full_aggr_signature, file = paste0(sig_data_path,"full/","full_aggregated_signature.tsv"))
    print(paste0("The full aggregated signature was saved at ",output_dir))

    # convert the full aggregated signature to matrix with row names: GeneID and col: aggregated_GeneScores for connectivity metric calculation
    row.names(full_aggr_signature) <- full_aggr_signature$GeneID
    full_aggr_signature$GeneID <- NULL
    full_aggr_sig_matrix <- as.matrix(full_aggr_signature)
    # get drug candidates using the specified Correlation-based method
    final_aggr_res <- compute_Cor_based_scores(full_aggr_sig_matrix, score_method, db_path)
}
write_tsv(final_aggr_res, file=paste0(output_dir,"/",score_method,"_",platform,"_aggregated_signature.tsv"))
print(paste0("Aggregated-signature drug results saved for ",sig_metadata_path," at ",output_dir,"/",score_method,"_",platform,"_aggregated_signature.tsv"))

## --------- individual signatures ----------- ##

# read in the input signature metadata
data_to_run <- read.delim(sig_metadata_path, sep="\t")

# iterate through each signature info
# compare the disease signature against drug profiles using the method of choice
for(i in 1:nrow(data_to_run)){
    print(paste0("iteration ", i))
    print(paste0("Method: ",score_method))
    # get info from the metadata file
    accession_no <- data_to_run$accession_no[i]
    file_name <- data_to_run$file_name[i]
    print(paste0("Quantifying drug candidates for ", accession_no," ",file_name))

    if(score_method %in% c("LINCS","CMAP")){
        # get signature paths
        up_sig_path <- paste0(sig_data_path,"up/",accession_no,"_",file_name,"_up.tsv")
        dn_sig_path <- paste0(sig_data_path,"dn/",accession_no,"_",file_name,"_dn.tsv")

        # get disease signatures
        if(file.exists(up_sig_path) && file.exists(dn_sig_path)){ # check if the signatures exist
            # up signature
            up_genes <- get_updn_signature(up_sig_path)
            # dn signature
            dn_genes <- get_updn_signature(dn_sig_path)
        }else{
            next
        }

    }else if(score_method %in% c("Cor_spearman","Cor_pearson")){
        full_sig_path <- paste0(sig_data_path,"full/",accession_no,"_",file_name,"_full.tsv")
        if(file.exists(full_sig_path)){
            full_signature_matrix <- get_full_signature(full_sig_path)
        }else{
            next
        }

    }
    # query the input gene signatures against all drug signatures in the selected database
    # finalize the drug result dataframe: add metadata info, filter out non-FDA-approved drugs, and remove rows with NAs
    if(score_method == "CMAP"){
        final_res <- compute_CMap1_scores(up_genes, dn_genes, db_path)
    }else if (score_method == "LINCS"){
        final_res <- compute_CMap2lincs_scores(up_genes, dn_genes, db_path)
    }else if(score_method %in% c("Cor_spearman","Cor_pearson")){
        final_res <- compute_Cor_based_scores(full_signature_matrix, score_method, db_path)
    }

    write_tsv(final_res, file=paste0(output_dir,"/",score_method,"_",accession_no,"_",file_name,".tsv"))
    print(paste0("Individual signature results saved for ", score_method,"_",accession_no,"_",file_name))
}

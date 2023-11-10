import pandas as pd
import numpy as np
import scipy.stats as ss
import rbo

def rm_na(arr1, arr2):
    arr1, arr2 = np.array(arr1), np.array(arr2)
    mask = ~(np.isnan(arr1) | np.isnan(arr2))
    return arr1[mask], arr2[mask]

def calculate_correlations(subset_drug_data, dis_data):
    correlation_matrix = pd.DataFrame(index=subset_drug_data.columns, columns=dis_data.columns)
    for drug_sample in subset_drug_data.columns:
        sorted_drug_geneids = subset_drug_data.sort_values(by=drug_sample, ascending=False).index
        for dis_sample in dis_data.columns:
            entry = np.array([np.nan, np.nan, np.nan, np.nan], dtype=float)
            # pearson
            corr = ss.pearsonr(*rm_na(subset_drug_data[drug_sample], dis_data[dis_sample]))[0]
            entry[0] = corr
            # spearman
            corr = ss.spearmanr(*rm_na(subset_drug_data[drug_sample], dis_data[dis_sample]))[0]
            entry[1] = corr
            # rbo
            sorted_dis_geneids = dis_data.sort_values(by=dis_sample, ascending=False).index
            corr = rbo.RankingSimilarity(*rm_na(sorted_drug_geneids, sorted_dis_geneids)).rbo()
            entry[2] = corr
            # rbo p=.99
            corr = rbo.RankingSimilarity(*rm_na(sorted_drug_geneids, sorted_dis_geneids)).rbo(p=.99)
            entry[3] = corr
            correlation_matrix.loc[drug_sample, dis_sample] = entry
    return correlation_matrix

def normalize_matrix(matrix):
    matrix_np = matrix.to_numpy(copy=True, dtype=float)

    row_means = np.mean(matrix_np, axis=1)
    row_stds = np.std(matrix_np, axis=1)
    col_means = np.mean(matrix_np, axis=0)
    col_stds = np.std(matrix_np, axis=0)

    zr = (matrix_np - row_means[:, np.newaxis]) / row_stds[:, np.newaxis]
    zc = (matrix_np - col_means) / col_stds

    result = (zr + zc) / np.sqrt(2)

    normalized_matrix = matrix.copy()
    normalized_matrix.loc[:, :] = result

    return normalized_matrix.astype(float)
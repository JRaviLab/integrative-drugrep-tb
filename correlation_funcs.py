import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as ss
from rbo import RankingSimilarity as rbo_rank
import gzip
from sklearn.linear_model import Lasso

def rm_na(arr1, arr2):
    arr1, arr2 = np.array(arr1), np.array(arr2)
    mask = ~(np.isnan(arr1) | np.isnan(arr2))
    return arr1[mask], arr2[mask]

def get_correlations(drug_data, dis_sample, correlation_type='pearson'):
    """
    Calculate correlations between a disease sample and each drug sample.

    Parameters
    ----------
    drug_data : pandas.DataFrame
        Expression data of drug samples
    dis_sample : pandas.Series
        Expression data of a disease sample
    correlation_type : str or array_like, optional
        Type of correlation to calculate. Default is 'pearson'.
        Options are 'pearson', 'spearman', 'rbo', 'lasso', 'all'.

    Returns
    -------
    output : pandas.DataFrame
        Correlations between the disease sample and each drug sample. Each index 
        is a drug sample and each column is a correlation type.
    """
    
    if correlation_type == 'all':
        correlation_type = ['pearson', 'spearman', 'rbo', 'lasso']
    correlation_type = np.array(correlation_type)

    output = pd.DataFrame(index=drug_data.columns, columns=correlation_type)

    if 'pearson' in correlation_type:
        for drug_sample in drug_data.columns:
            output.loc[drug_sample, 'pearson'] = \
                ss.pearsonr(*rm_na(dis_sample, drug_data[drug_sample]))[0]
    if 'spearman' in correlation_type:
        for drug_sample in drug_data.columns:
            output[drug_sample, 'spearman'] = \
                ss.spearmanr(*rm_na(dis_sample, drug_data[drug_sample]))[0]
    if 'lasso' in correlation_type:
        dis_sample.dropna(inplace=True)
        drug_data_nona = drug_data.loc[dis_sample.index]
        lasso = Lasso(alpha=0.1)
        lasso.fit(drug_data_nona.values, dis_sample.values)
        output[drug_data_nona.index, 'lasso'] = lasso.coef_
    if 'rbo' in correlation_type:
        dis_sample.dropna(inplace=True)
        sorted_dis_geneids = dis_sample.sort_values(ascending=False).index
        drug_data = drug_data[np.isin(drug_data.index, sorted_dis_geneids)]
        for drug_sample in drug_data.columns:
            sorted_drug_geneids = drug_data[drug_sample].sort_values(ascending=False).index
            output[drug_sample, 'rbo'] = \
                rbo_rank(*rm_na(sorted_dis_geneids, sorted_drug_geneids)).rbo()
            
    return output

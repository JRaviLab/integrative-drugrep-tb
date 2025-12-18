import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as ss
from rbo import RankingSimilarity as rbo_rank
import gzip
from sklearn.linear_model import Lasso
import h5py

def rm_na(arr1, arr2):
    arr1, arr2 = np.array(arr1), np.array(arr2)
    mask = ~(np.isnan(arr1) | np.isnan(arr2))
    return arr1[mask], arr2[mask]

def get_correlations(drug_data, dis_sample, correlation_type='pearson'):
    """
    Calculate correlations between a disease sample and each drug sample.
    Expects no missing values in the drug data or disease sample.

    Parameters
    ----------
    drug_data : pandas.DataFrame
        Expression data of drug samples
    dis_sample : pandas.Series or pandas.DataFrame
        Expression data of disease sample(s)
    correlation_type : str or array_like, optional
        Type of correlation to calculate. Default is 'pearson'.
        Options are 'pearson', 'spearman', 'rbo', 'lasso', 'all'.

    Returns
    -------
    output : pandas.DataFrame
        Correlations between the disease sample and each drug sample. Each index 
        is a drug sample and each column is a correlation type.
    """

    correlation_type = np.array(correlation_type)
    if 'all' in correlation_type :
        correlation_type = np.array(['pearson', 'spearman', 'rbo', 'lasso'])


    if isinstance(dis_sample, pd.DataFrame):
        # array_3d[dis_sample, drug_sample, correlation_type]
        array_3d = np.array([
            get_correlations(drug_data, dis_sample[col], correlation_type).values
            for col in dis_sample.columns
        ], dtype=float)
        return array_3d

    else:
        output = pd.DataFrame(index=drug_data.columns, columns=correlation_type)

        if 'pearson' in correlation_type:
            for drug_sample in drug_data.columns:
                output.loc[drug_sample, 'pearson'] = \
                    ss.pearsonr(dis_sample, drug_data[drug_sample])[0]
        if 'spearman' in correlation_type:
            for drug_sample in drug_data.columns:
                output.loc[drug_sample, 'spearman'] = \
                    ss.spearmanr(dis_sample, drug_data[drug_sample])[0]
        if 'rbo' in correlation_type:
            sorted_dis_geneids = dis_sample.sort_values(ascending=False).index.values
            drug_data = drug_data[np.isin(drug_data.index, sorted_dis_geneids)]
            for drug_sample in drug_data.columns:
                sorted_drug_geneids = drug_data[drug_sample].sort_values(ascending=False).index.values
                output.loc[drug_sample, 'rbo'] = \
                    rbo_rank(sorted_dis_geneids, sorted_drug_geneids).rbo(p=.99) 
        if 'lasso' in correlation_type:
            lasso = Lasso(alpha=.1, max_iter=10000)
            lasso.fit(drug_data.values, dis_sample.values)
            output.loc[:, 'lasso'] = lasso.coef_
     
        return output

def get_cell_line_vector(correlation_df):
    """
    Calculate the cell line vector of a correlation matrix.
    
    Parameters
    ----------
    correlation_df : pandas.DataFrame
        Correlations for each drug sample. Each index is a drug sample and each 
        column is a correlation type.
        
    Returns
    -------
    output : pandas.Series
        Cell line vector of the correlation matrix.
    """
    cell_lines = pd.Series(index=correlation_df.index, dtype=str)
    for drug_sample in correlation_df.index:
        cell_lines[drug_sample] = drug_sample.split('_')[1]
    unique_cell_lines = np.unique(cell_lines)

    output = pd.Series(index=unique_cell_lines)
    for drug_sample in correlation_df.index:
        output[drug_sample] = np.mean(correlation_df.loc[drug_sample])
    return output

def cell_line_ranks(pairwise_df):
    """
    Calculate a df of ranked cell lines for each disease sample from pairwise
    correlations.
    
    Parameters
    ----------
    pairwise_df : pandas.DataFrame
        Correlations for each sample pair. Each index is a drug sample and each
        column is a disease sample.
    
    Returns
    -------
    output : pandas.DataFrame
        Columns are disease samples and rows are cell lines in order of average
    """
    cell_lines = pd.Series(index=pairwise_df.index, dtype=str)
    for drug_sample in pairwise_df.index:
        cell_lines[drug_sample] = drug_sample.split('_')[1]
    unique_cell_lines = np.unique(cell_lines)

    output = {}
    for dis_sample in pairwise_df.columns:
        cell_line_means = pd.Series(index=unique_cell_lines)
        for cell_line in unique_cell_lines:
            cell_line_samples = cell_lines[cell_lines == cell_line].index
            cell_line_means[cell_line] = pairwise_df.loc[cell_line_samples, dis_sample].mean()
        cell_line_means.sort_values(ascending=False, inplace=True)
        output[dis_sample] = cell_line_means.index.values
    return pd.DataFrame(output.values(), index=output.keys()).T

def avg_rank(dfs):
    """
    Calculate one dataframe by comparing average of input dataframes with 
    ranked lists as columns.
    
    Parameters
    ----------
    dfs : array-like of pandas.DataFrame
        List of dataframes to average.
        
    Returns
    -------
    output : pandas.DataFrame
        Dataframe with averaged ranked lists.
    """
    
    if len(dfs) < 2:
        raise ValueError('Must input at least 2 dataframes')
    
    if not all(dfs[0].shape == df.shape for df in dfs):
        raise ValueError('Dataframes must all have the same shape')
    
    output = pd.DataFrame(index=dfs[0].index, columns=dfs[0].columns)
    for col in dfs[0].columns:
        ranks = pd.Series(0, index=dfs[0].iloc[:, 0])
        for row in dfs[0].index:
            for df in dfs:
                ranks[df.loc[row, col]] += row
        output[col] = ranks.sort_values(ascending=True).index
    return output
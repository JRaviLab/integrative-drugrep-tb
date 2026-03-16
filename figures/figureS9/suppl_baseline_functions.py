from collections import defaultdict
import os
import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from rbo import RankingSimilarity as rbo_rank

mpl.rcParams['font.family'] = 'Arial'

# ///// Constants 
ALL_CORRELATION_TYPES = ['pearson', 'spearman', 'rbo', 'lasso']
RBO_P                 = 0.99
LASSO_ALPHA           = 0.1
LASSO_MAX_ITER        = 10_000
STOUFFER_MIN_CORR     = 0.05

# Label shortening maps
DIS_TISSUE_NAME_MAP = {
    'primary human lymphatic endothelial cell': 'lymphatic endothelial',
    'PBMCs isolated from Healthy individuals' : 'PBMC',
    'Purified cord blood CD34+ cells'         : 'cord blood CD34+',
    'monocyte-derived macrophages'            : 'macrophages',
}
DRUG_TISSUE_NAME_MAP = {
    'haematopoietic and lymphoid tissue': 'haematopoietic\n& lymphoid',
    'central nervous system'            : 'central nervous\nsystem',
}

# ///// Data loading & preprocessing

def load_drug_data(path: str, time_filter: str = '24H') -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load L1000 Level 3 drug baseline data, filter to a single time point,
    and return the expression matrix alongside a parsed sample metadata table.

    Returns
    -------
    drug_data : pd.DataFrame  genes × samples (filtered)
    split_drug : pd.DataFrame  parsed sample metadata
    """
    raw = pd.read_csv(path, sep='\t', index_col='rid')

    split = pd.DataFrame(
        [s.split('_', 4) for s in raw.columns],
        columns=['exp', 'cell', 'time', 'x', 'end'],
    )
    split = split[split['time'] == time_filter].sort_values(['exp', 'cell', 'end'])
    split.reset_index(drop=True, inplace=True)

    sample_ids = (
        split['exp'] + '_' + split['cell'] + '_' +
        split['time'] + '_' + split['x'] + '_' + split['end']
    )
    drug_data = raw[sample_ids]

    return drug_data, split

def average_replicates(drug_data: pd.DataFrame, split_drug: pd.DataFrame) -> pd.DataFrame:
    """
    Average expression across technical replicates that share the same
    (exp, cell, time, end) key, collapsing the dose ('x') dimension.
    """
    def sample_key(row, include_x=False):
        cols = ['exp', 'cell', 'time'] + (['x'] if include_x else []) + ['end']
        return '_'.join(str(row[col]) for col in cols)

    split_drug = split_drug.copy()
    split_drug['key']   = split_drug.apply(sample_key, axis=1)
    split_drug['key_x'] = split_drug.apply(lambda r: sample_key(r, include_x=True), axis=1)

    averaged = {}
    for key, grp in split_drug.groupby('key', sort=False):
        sample_ids = grp['key_x'].tolist()
        averaged[key] = drug_data[sample_ids].mean(axis=1)

    return pd.DataFrame(averaged)

def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Quantile-normalize all columns to the same distribution."""
    sorted_vals = np.sort(df.values, axis=0)
    rank_means  = sorted_vals.mean(axis=1)
    result = df.copy()
    for col in df.columns:
        ranks = np.searchsorted(np.sort(df[col]), df[col])
        result[col] = rank_means[ranks]
    return result

def stouffer_aggregate(drug_data: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse per-drug samples into one weighted vector per cell line using
    Stouffer-like weighting (pairwise inter-replicate correlations as weights).
    Falls back to simple mean when no replicates agree above the threshold.
    """
    cell_lines = sorted({col.split('_')[1] for col in drug_data.columns})
    aggregated  = pd.DataFrame(index=drug_data.index, columns=cell_lines, dtype=float)

    for cell_line in cell_lines:
        cols   = [c for c in drug_data.columns if f'_{cell_line}_' in c]
        subset = drug_data[cols].values          # genes × replicates

        if len(cols) == 1:
            aggregated[cell_line] = subset[:, 0]
            continue

        corr_matrix = np.corrcoef(subset.T)      # replicates × replicates
        n    = len(cols)
        mask = ~np.eye(n, dtype=bool)
        weights = corr_matrix[mask].reshape(n, n - 1).mean(axis=1)
        weights = np.where(weights < STOUFFER_MIN_CORR, 0.0, weights)

        weight_sum = weights.sum()
        aggregated[cell_line] = (
            subset.mean(axis=1) if weight_sum == 0
            else subset @ (weights / weight_sum)
        )

    return aggregated

def load_disease_data(
    rnaseq_dir: str,
    micro_dir: str,
    rna_meta_path: str,
    micro_meta_path: str,
    gene_info_path: str,
    landmark_genes: pd.Index,
    drug_reference_col: pd.Series,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """
    Load and preprocess RNAseq and Microarray disease control samples:
      - Filter to healthy controls
      - Subset to landmark genes
      - Quantile-normalize to match drug data distribution

    Returns
    -------
    dis_data_qn : pd.DataFrame  genes x samples, QN-matched to drug data
    tissue       : pd.Series    tissue label per sample
    cell         : pd.Series    cell label per sample
    """
    rna_meta = pd.read_csv(rna_meta_path, sep='\t')
    micro_meta = pd.read_csv(micro_meta_path, sep='\t')
    meta = pd.concat([rna_meta, micro_meta])
    controls = meta.loc[meta['CLASSIFICATION'] == 'healthy control without treatment', 'geo_accession'].values

    gene_info = pd.read_csv(gene_info_path, sep='\t', header=0, index_col=0)
    gene_info = gene_info[[g in landmark_genes for g in gene_info['GeneID'].values]]
    gene_map = gene_info.drop_duplicates('Ensembl').set_index('Ensembl')['GeneID']

    dis_data = pd.DataFrame({'GeneID': landmark_genes})
    for fname in os.listdir(rnaseq_dir):
        exp = pd.read_csv(
            os.path.join(rnaseq_dir, fname),
            sep='\t', low_memory=False, dtype={'GeneID': 'int32'},
        )
        exp.columns.values[0] = 'GeneID'
        exp['GeneID'] = exp['GeneID'].map(gene_map)
        exp = exp[exp['GeneID'].isin(landmark_genes)]
        exp = exp[exp.columns[exp.columns.isin(np.append(controls, 'GeneID'))]]
        dis_data = dis_data.merge(exp, how='left', on='GeneID')
    for fname in os.listdir(micro_dir):
        exp = pd.read_csv(
            os.path.join(micro_dir, fname),
            sep='\t', low_memory=False, dtype={'GeneID': 'int32'},
        )
        exp.columns.values[0] = 'GeneID'
        exp['GeneID'] = exp['GeneID'].map(gene_map)
        exp = exp[exp['GeneID'].isin(landmark_genes)]
        exp = exp[exp.columns[exp.columns.isin(np.append(controls, 'GeneID'))]]
        dis_data = dis_data.merge(exp, how='left', on='GeneID')

    dis_data = dis_data.set_index('GeneID').astype(float)
    dis_data.columns.name = 'Disease Samples'
    dis_data = dis_data.T.drop_duplicates().T
    dis_data.dropna(inplace=True)

    # Quantile-normalize disease samples to the drug data distribution
    target = np.sort(drug_reference_col.loc[dis_data.index].values)
    dis_qn = dis_data.copy()
    for col in dis_data.columns:
        ranks = np.searchsorted(np.sort(dis_data[col]), dis_data[col])
        dis_qn[col] = target[ranks]

    tissue = meta.set_index('geo_accession').loc[dis_data.columns, 'TISSUE'].fillna('Missing')
    cell = meta.set_index('geo_accession').loc[dis_data.columns, 'CELL'].fillna('Missing')

    return dis_qn.astype(float), tissue, cell

# ///// Correlation computation

def get_correlations(
    drug_data: pd.DataFrame,
    dis_sample: pd.Series | pd.DataFrame,
    correlation_type: str | list[str] = 'pearson',
) -> pd.DataFrame | np.ndarray:
    """
    Calculate correlations between disease sample(s) and each drug cell line.

    Parameters
    ----------
    drug_data        : genes × cell-lines
    dis_sample       : genes (Series) or genes × samples (DataFrame)
    correlation_type : one or more of 'pearson','spearman','rbo','lasso','all'

    Returns
    -------
    pd.DataFrame  shape (n_cell_lines, n_types)   — when dis_sample is a Series
    np.ndarray    shape (n_samples, n_cell_lines, n_types) — for a DataFrame
    """
    types = np.atleast_1d(correlation_type)
    if 'all' in types:
        types = np.array(ALL_CORRELATION_TYPES)

    if isinstance(dis_sample, pd.DataFrame):
        return np.array(
            [get_correlations(drug_data, dis_sample[col], types).values
             for col in dis_sample.columns],
            dtype=float,
        )

    output = pd.DataFrame(index=drug_data.columns, columns=types, dtype=float)

    if 'pearson' in types:
        output['pearson'] = [
            ss.pearsonr(dis_sample, drug_data[col])[0]
            for col in drug_data.columns
        ]
    if 'spearman' in types:
        output['spearman'] = [
            ss.spearmanr(dis_sample, drug_data[col])[0]
            for col in drug_data.columns
        ]
    if 'rbo' in types:
        dis_ranked = dis_sample.sort_values(ascending=False).index.values
        rbo_data   = drug_data.loc[drug_data.index.isin(dis_ranked)]
        output['rbo'] = [
            rbo_rank(dis_ranked,
                     rbo_data[col].sort_values(ascending=False).index.values
                     ).rbo(p=RBO_P)
            for col in rbo_data.columns
        ]
    if 'lasso' in types:
        lasso = Lasso(alpha=LASSO_ALPHA, max_iter=LASSO_MAX_ITER)
        lasso.fit(drug_data.values, dis_sample.values)
        output['lasso'] = lasso.coef_

    return output

def compute_correlation_matrix(
    drug_data: pd.DataFrame,
    dis_data: pd.DataFrame,
    correlation_types: list[str],
) -> dict[str, pd.DataFrame]:
    """
    Compute a correlation matrix between all dis_data columns and all
    drug_data columns for each requested correlation type.

    This is the single entry point used for BOTH drug×disease and drug×drug —
    ensuring identical ordering and aggregation logic in both cases.

    Returns
    -------
    dict  corr_type → pd.DataFrame  shape (n_drug_cells, n_dis_samples)
    """
    # 3-D array: (n_dis_samples, n_drug_cells, n_types)
    arr = get_correlations(drug_data, dis_data, correlation_types).astype('float32')
    results = {}
    for i, ctype in enumerate(correlation_types):
        # rows = drug cell lines, cols = disease samples
        results[ctype] = pd.DataFrame(
            arr[:, :, i].T,
            index=drug_data.columns,
            columns=dis_data.columns,
        )
    return results

# ///// Tissue-level aggregation

def build_tissue_groups(
    series: pd.Series,
    min_count: int = 0,
) -> dict[str, list]:
    """
    Build a tissue → [sample/cell-line ids] mapping from a Series whose
    index is sample/cell ids and values are tissue labels.
    Optionally filter to tissues with at least `min_count` members.
    """
    groups: dict[str, list] = defaultdict(list)
    counts = series.value_counts()
    for item, tissue in series.items():
        if counts[tissue] > min_count:
            groups[tissue].append(item)
    return dict(groups)

def build_tissue_matrix(
    corr_matrix: pd.DataFrame,
    row_groups: dict[str, list],
    col_groups: dict[str, list],
    mask_diagonal: bool = False,
) -> pd.DataFrame:
    """
    Aggregate a sample-level correlation matrix into a tissue × tissue matrix
    of median correlations.

    Parameters
    ----------
    corr_matrix    : rows = drug cell lines, cols = disease samples (or vice versa)
    row_groups     : tissue → [row ids]
    col_groups     : tissue → [col ids]
    mask_diagonal  : if True, NaN-out diagonal cells (for square drug×drug matrix)
    """
    row_tissues = sorted(row_groups.keys())
    col_tissues = sorted(col_groups.keys())
    result = pd.DataFrame(index=row_tissues, columns=col_tissues, dtype=float)

    for rt in row_tissues:
        for ct in col_tissues:
            sub = corr_matrix.loc[
                corr_matrix.index.intersection(row_groups[rt]),
                corr_matrix.columns.intersection(col_groups[ct]),
            ].copy()

            if mask_diagonal and rt == ct:
                np.fill_diagonal(sub.values, np.nan)

            result.loc[rt, ct] = (
                np.nanmedian(sub.values)
                if np.isfinite(sub.values).any()
                else np.nan
            )

    return result

# ///// Visualization

def safe_zscore(matrix: np.ndarray, axis: int) -> np.ndarray:
    mean = np.nanmean(matrix, axis=axis, keepdims=True)
    std  = np.nanstd(matrix,  axis=axis, ddof=1, keepdims=True)
    return (matrix - mean) / std

def combined_zscore(matrix: np.ndarray) -> np.ndarray:
    """Row + column z-scores averaged (scaled by √2)."""
    return (safe_zscore(matrix, axis=1) + safe_zscore(matrix, axis=0)) / np.sqrt(2)

def make_axis_labels(
    tissues: list[str],
    groups: dict[str, list],
    name_map: dict[str, str],
    count_source: dict[str, int] | None = None,
) -> list[str]:
    """
    Build axis tick labels with pretty names and member counts.
    `count_source` overrides the group length when sample counts differ
    from cell-line counts (e.g., total drug samples vs. cell lines).
    """
    labels = []
    for t in tissues:
        pretty = name_map.get(t, t)
        count  = count_source[t] if count_source else len(groups[t])
        labels.append(f"{pretty} ({count})")
    return labels

def plot_heatmap(
    matrix: np.ndarray,
    xlabels: list[str],
    ylabels: list[str],
    title: str,
    cbar_label: str,
    xlabel: str = 'Tissue Type',
    ylabel: str = 'Tissue Type',
    cmap: str = 'Blues',
    annotate: bool = True,
    center_zero: bool = False,
) -> None:
    """
    Plot an annotated heatmap. NaN cells are shown in black.
    Text colour switches to white above the 75th percentile (or 75 % of absmax).
    """
    fig, ax = plt.subplots(figsize=(9, 8))

    if center_zero:
        absmax = np.nanmax(np.abs(matrix)) - 0.5
        vmin, vmax = -absmax, absmax
    else:
        absmax = None
        vmin, vmax = np.nanmin(matrix), np.nanmax(matrix)

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(color='black')

    cax  = ax.imshow(matrix, cmap=cmap_obj, aspect='auto', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cax, ax=ax, orientation='vertical', fraction=0.04, pad=0.04)
    cbar.set_label(cbar_label, fontsize=11)

    if annotate:
        threshold = absmax * 0.75 if center_zero else np.nanpercentile(matrix, 75)
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix[i, j]
                if not np.isnan(val):
                    contrast = abs(val) < threshold if center_zero else val < threshold
                    color = 'black' if contrast else 'white'
                    ax.text(j, i, f"{val:.2f}", ha='center', va='center',
                            fontsize=8, color=color)

    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=15, pad=20)

    plt.savefig(f'{title}.png', dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()

def plot_tissue_heatmaps(
    tissue_matrix: pd.DataFrame,
    xlabels: list[str],
    ylabels: list[str],
    title_prefix: str,
    xlabel: str,
    ylabel: str,
) -> None:
    """Plot both the raw median and z-scored heatmaps for a tissue matrix."""
    # Drop first row/col only for drug×drug (square) to remove the
    # lowest-count tissue; for rectangular matrices pass the full frame.
    mat = tissue_matrix.values.astype(float)

    plot_heatmap(
        mat, xlabels, ylabels,
        f'Median Pearson Correlation — {title_prefix}',
        'Pearson Corr.',
        xlabel=xlabel, ylabel=ylabel,
    )
    plot_heatmap(
        combined_zscore(mat), xlabels, ylabels,
        f'Z-Score Median Pearson Correlation — {title_prefix}',
        'Z-Score of Pearson Corr.',
        xlabel=xlabel, ylabel=ylabel,
        cmap='RdBu_r', center_zero=True,
    )

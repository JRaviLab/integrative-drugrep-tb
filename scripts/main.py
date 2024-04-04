import data_funcs
import correlation_funcs as cf
import pandas as pd
import time
import numpy as np
import shelve
import argparse
from datetime import date


def main(p=False, s=False, r=False, l=False, group=None):

    if not p and not s and not r and not l:
        print("No correlation type selected. Please select at least one of the following: pearson, spearman, rbo, lasso")
        return
    
    if group not in ['avg_cell_line', 'stouffer', 'none']:
        print("Invalid group option. Please select one of the following: avg_cell_line, stouffer, none")
        return
    
    start = time.time()

    # Dataframes: [gene, sample]
    # drug_data, landmark_genes = data_funcs.read_drug_data(quantile_normalize=True)
    qn_drug_lv3data = shelve.open("data/expression/drug_data/qn_drug_lv3data") 
    qndrug_target_vec_ = qn_drug_lv3data['qndrug_target_vec']
    drug_data_qn_ = qn_drug_lv3data['drug_data_qn']
    drug_data = drug_data_qn_
    landmark_genes = drug_data.index.values
    rnaseq_dis_data, marray_dis_data = data_funcs.read_dis_data(landmark_genes)

    if group == 'avg_cell_line':
        # group drug samples by cell line -> drug_data_cell
        cell_lines = list(set([x.split('_')[1] for x in drug_data.columns]))
        drug_data_cell = pd.DataFrame(index=drug_data.index, columns=cell_lines)
        for cell_line in cell_lines:
            cell_line_columns = [x for x in drug_data.columns if cell_line in x]
            drug_data_cell[cell_line] = drug_data[cell_line_columns].mean(axis=1)
        drug_data = drug_data_cell
    elif group == 'stouffer':
        # group drug samples by cell line -> drug_data_cell
        cell_lines = list(set([x.split('_')[1] for x in drug_data.columns]))
        drug_data_cell = pd.DataFrame(index=drug_data.index, columns=cell_lines)
        for cell_line in cell_lines:
            cell_line_columns = [x for x in drug_data.columns if cell_line in x]
            weights = []
            for col in cell_line_columns:
                correlations = []
                for col2 in cell_line_columns:
                    if col != col2:
                        correlations.append(np.corrcoef(drug_data[col], drug_data[col2])[0][1])
                avg_corr = np.mean(correlations)
                avg_corr = 0 if avg_corr<0.05 else avg_corr
                weights.append(np.mean(correlations))
            weights = np.array(weights)
            weights = weights/np.sum(weights)
            drug_data_cell[cell_line] = np.sum(drug_data[cell_line_columns].values * weights, axis=1)

        drug_data = drug_data_cell


    def qnorm_dis_data(dis_data, target_vec):
        dis_data_qn = dis_data.copy()
        target_vec = np.sort(np.array(target_vec))
        for col in dis_data_qn.columns:
            t = np.searchsorted(np.sort(dis_data[col]), dis_data[col])
            dis_data_qn.loc[:,col] = [target_vec[i] for i in t]
        return dis_data_qn.astype('float64')

    rnaseq_dis_data.dropna(inplace=True)
    rnaseq_genes = rnaseq_dis_data.index
    rnaseq_dis_data_qn = qnorm_dis_data(rnaseq_dis_data, qndrug_target_vec_.loc[rnaseq_genes])
    marray_dis_data.dropna(inplace=True)
    marray_genes = marray_dis_data.index
    marray_dis_data_qn = qnorm_dis_data(marray_dis_data, qndrug_target_vec_)

    types = ['pearson', 'spearman', 'rbo', 'lasso']
    mask = [p, s, r, l]
    types = np.array(types)[mask]

    def save_ranks(drug_data, dis_data, filename):
        # 3D array: [drug, dis, corr]
        arr_3d = cf.get_correlations(drug_data, dis_data, 
                                    correlation_type=types).astype('float32')
        
        # List of dataframes: [drug, dis]
        corr_matrices = [pd.DataFrame(arr_3d[:,:,i].T, index=drug_data.columns, 
                                    columns=dis_data.columns) 
                        for i in range(arr_3d.shape[2])]
        
        for i, mat in enumerate(corr_matrices):
            mat.to_csv(f'{filename}{types[i][0]}_{str(date.today())[5:]}.csv')

    drug_data_rna_subset = drug_data.loc[rnaseq_genes]
    drug_data_marray_subset = drug_data.loc[marray_genes]

    save_ranks(drug_data_rna_subset, rnaseq_dis_data_qn, f'results/rna_ranks_')
    save_ranks(drug_data_marray_subset, marray_dis_data_qn, f'results/marray_ranks_')

    print('Time (s):', time.time()-start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Description of your script")
    parser.add_argument('-p', '--pearson', action="store_true", help="Caclulate pearson correlation")
    parser.add_argument('-s', '--spearman', action="store_true", help="Caclulate spearman correlation")
    parser.add_argument('-r', '--rbo', action="store_true", help="Caclulate rbo correlation")
    parser.add_argument('-l', '--lasso', action="store_true", help="Caclulate lasso betas")
    parser.add_argument('-g', '--group', type=str, default=None, help="How to group drug samples. Options: avg_cell_line, stouffer, none. Default: none")
    
    args = parser.parse_args()
    
    main(p=args.pearson, s=args.spearman, r=args.rbo, l=args.lasso, group=args.group)
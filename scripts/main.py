import data_funcs
import correlation_funcs as cf
import pandas as pd
import time
import numpy as np
import shelve
import argparse


def main(p=False, s=False, r=False, l=False):
    start = time.time()

    # Dataframes: [gene, sample]
    # drug_data, landmark_genes = data_funcs.read_drug_data(quantile_normalize=True)
    qn_drug_lv3data = shelve.open("data/expression/drug_data/qn_drug_lv3data") 
    qndrug_target_vec_ = qn_drug_lv3data['qndrug_target_vec']
    drug_data_qn_ = qn_drug_lv3data['drug_data_qn']
    drug_data = drug_data_qn_
    landmark_genes = drug_data.index.values
    rnaseq_dis_data, marray_dis_data = data_funcs.read_dis_data(landmark_genes)

    # group drug samples by cell line -> drug_data_cell
    cell_lines = list(set([x.split('_')[1] for x in drug_data.columns]))
    drug_data_cell = pd.DataFrame(index=drug_data.index, columns=cell_lines)
    for cell_line in cell_lines:
        cell_line_columns = [x for x in drug_data.columns if cell_line in x]
        drug_data_cell[cell_line] = drug_data[cell_line_columns].mean(axis=1)
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

        # if len(corr_matrices) == 1:
        #     output = corr_matrices[0]
        # else:
        #     # Dataframe: [drug, dis]
        #     output = cf.avg_rank([cf.cell_line_ranks(mat) for mat in corr_matrices])
        
        for i, mat in enumerate(corr_matrices):
            mat.to_csv(f'{filename}{types[i][0]}.csv')

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
    
    args = parser.parse_args()
    
    main(p=args.pearson, s=args.spearman, r=args.rbo, l=args.lasso)
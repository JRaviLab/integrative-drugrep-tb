import landon_scripts.data_funcs_old as data_funcs_old
import correlation_funcs as cf
import pandas as pd
import time

start = time.time()

# Dataframes: [gene, sample]
drug_data, landmark_genes = data_funcs_old.read_drug_data()
dis_data = data_funcs_old.read_dis_data(landmark_genes)

# 3D array: [drug, dis, corr]
arr_3d = cf.get_correlations(drug_data, dis_data, 
                             correlation_type=['pearson', 'spearman', 'rbo']).astype('float32')

# List of dataframes: [drug, dis]
corr_matrices = [pd.DataFrame(arr_3d[:,:,i].T, index=drug_data.columns, 
                              columns=dis_data.columns) 
                 for i in range(arr_3d.shape[2])]

# Dataframe: [drug, dis]
output = cf.avg_rank([cf.cell_line_ranks(mat) for mat in corr_matrices])

output.to_csv('all_combined_cell_line_ranks.csv', index=False)

print('Time (s):', time.time()-start)
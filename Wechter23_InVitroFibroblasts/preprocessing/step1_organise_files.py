
import scanpy as sc
import pandas as pd

data_directory = '../data/Wechter2023/''

import os

srr_runs = pd.read_csv(data_directory + 'SRR_Acc_List_static.txt', header=None)[0].tolist()

srr_run_metadata = pd.read_csv(data_directory + 'Wechter2023/SraRunTable.txt')
srr_run_metadata = srr_run_metadata[srr_run_metadata['Run'].isin(srr_runs)]

for srr_run in srr_runs:
    
    filepath = data_directory

    foldername = 'out_' + srr_run + '_v3'
    srr_run_meta = srr_run_metadata[srr_run_metadata['Run'] == srr_run]
    
    adata_filepath = os.path.join(filepath, foldername + '/counts_unfiltered/')
    adata_sample = sc.read_h5ad(adata_filepath + 'adata.h5ad')
    adata_gene_names = pd.read_csv(adata_filepath + 'cells_x_genes.genes.names.txt', header=None)
    adata_sample.var['gene_id'] = pd.Series(adata_sample.var_names, dtype='str').values
    adata_sample.var_names = pd.Index(adata_gene_names[0].values)

    # We'll define the counts as mature + ambiguous, i.e., exonic reads
    adata_sample.X = adata_sample.layers['ambiguous'] + adata_sample.layers['mature'] + adata_sample.layers['nascent']
    sc.pp.filter_genes(adata_sample, min_cells=1)
    sc.pp.filter_cells(adata_sample, min_genes=1)
    
    for col in srr_run_meta.columns:
        adata_sample.obs[col] = srr_run_meta[col].values[0]
    
    print(adata_sample)
    adata_sample.write(adata_filepath + 'adata_filtered.h5ad', compression='gzip')
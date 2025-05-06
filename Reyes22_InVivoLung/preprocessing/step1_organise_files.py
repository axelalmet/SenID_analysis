
import scanpy as sc
import pandas as pd

data_directory = '../data/Reyes2022/'

import os

srr_runs_gfppos = pd.read_csv(data_directory + 'SRR_Acc_List_Old_GFPPos.txt', header=None)[0].tolist()
srr_runs_gfpneg = pd.read_csv(data_directory + 'SRR_Acc_List_Old_GFPNeg.txt', header=None)[0].tolist()

srr_runs = srr_runs_gfppos + srr_runs_gfpneg

for srr_run in srr_runs:
    
    filepath = os.path.join(data_directory + 'Old_GFPPos/')
    if srr_run in srr_runs_gfpneg:
        filepath = os.path.join(data_directory + 'Old_GFPNeg/')

    foldername = 'out_' + srr_run + '_v3'
    
    adata_filepath = os.path.join(filepath, foldername + '/counts_unfiltered/')
    adata_sample = sc.read_h5ad(adata_filepath + 'adata.h5ad')
    adata_gene_names = pd.read_csv(adata_filepath + 'cells_x_genes.genes.names.txt', header=None)
    adata_sample.var['gene_id'] = pd.Series(adata_sample.var_names, dtype='str').values
    adata_sample.var_names = pd.Index(adata_gene_names[0].values)

    adata_sample.obs['Age'] = 'Old'

    gfp_status = 'Pos'

    if srr_run in srr_runs_gfpneg:
        gfp_status = 'Neg'

    adata_sample.obs['Age_GFP'] = 'Old_GFP' + gfp_status
    adata_sample.obs['Run'] = srr_run
    
    adata_sample.X = adata_sample.layers['ambiguous'] + adata_sample.layers['mature'] + adata_sample.layers['nascent']
    sc.pp.filter_genes(adata_sample, min_cells=1)
    sc.pp.filter_cells(adata_sample, min_genes=1)
        
    adata_sample.write(adata_filepath + 'adata_filtered.h5ad', compression='gzip')

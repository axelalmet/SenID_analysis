import scanpy as sc
from scipy.sparse import csr_matrix

samples = [f'P{i}' for i in range(1, 5)]

data_directory = '../data/Park2024/'

for sample in samples:
    adata = sc.read_h5ad(data_directory + sample + '/park24_' + sample + '_spe.h5ad')
    adata.var_names_make_unique()
    adata.X = csr_matrix(adata.X)
    adata.layers['counts'] = csr_matrix(adata.layers['counts'])
    adata.obsm['spatial'] = adata.obs[['array_col', 'array_row']].values
    
    adata.write(data_directory + sample + '/park24_' + sample + '_spe.h5ad', compression='gzip')


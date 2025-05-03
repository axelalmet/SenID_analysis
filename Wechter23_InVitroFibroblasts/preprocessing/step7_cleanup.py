
import scanpy as sc
from scipy.sparse import csr_matrix

data_directory = '../data/Wechter2023/'

adata = sc.read_h5ad(data_directory + 'wechter23_merged_unspliced.h5ad')

adata.obs_names_make_unique()
adata.var_names_make_unique()

# Make sure they're CSR sparse matrices (previous versions of zellkonverter set them to be CSC)
adata.X = csr_matrix(adata.X)
for layer in adata.layers:
    adata.layers[layer] = csr_matrix(adata.layers[layer])

adata.obs['Condition'] = 'CTRL'
adata.obs['Condition'][adata.obs['treatment'].str.contains('Replicative')] = 'RS'
adata.obs['Condition'][adata.obs['treatment'].str.contains('Etoposide')] = 'ETO'
adata.obs['Condition'][adata.obs['treatment'].str.contains('IR-induced')] = 'IR'

adata.write(data_directory + 'wechter23_merged_unspliced.h5ad', compression='gzip')  


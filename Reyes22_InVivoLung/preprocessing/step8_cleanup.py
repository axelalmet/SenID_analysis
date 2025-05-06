
import scanpy as sc
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (5, 4)
sc.settings.figdir = '/Users/axelalmet/Documents/Senescence/Figures/'
data_directory = '/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Reyes2022/'

adata = sc.read_h5ad(data_directory + 'reyes22_old_unspliced.h5ad')

adata.obs['GFP'] = 'Pos'
adata.obs['GFP'][adata.obs['Age_GFP'].str.contains('Neg')] = 'Neg'

adata.obs_names_make_unique()
adata.var_names_make_unique()

adata.X = csr_matrix(adata.X)
for layer in adata.layers:
    adata.layers[layer] = csr_matrix(adata.layers[layer])

adata.layers['spliced'] = adata.layers['mature'] + adata.layers['ambigous']
adata.layers['unspliced'] = adata.layers['nascent']

del adata.layers['mature'], adata.layers['ambigous'], adata.layers['nascent']

adata.write(data_directory + 'reyes22_old_unspliced.h5ad', compression='gzip')  

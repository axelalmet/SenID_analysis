import scvi
import torch
import scanpy as sc
import numpy as np

torch.set_float32_matmul_precision('high')

data_directory = '../data/Reyes2022/'

adata_ref = sc.read_h5ad(data_directory + 'tabula_muris_Lung_old.h5ad')
adata_query = sc.read_h5ad(data_directory + 'reyes22_old_unspliced.h5ad')

# Set the counts layer to be mature + ambiguous
adata_ref.layers['counts'] = adata_ref.X.copy()
adata_query.layers['counts'] = adata_query.layers['mature'] + adata_query.layers['ambiguous']

# Set which one is the reference and query
adata_ref.obs['label'] = 'reference'
adata_query.obs['label'] = 'query'

adata_query.obs['cell_ontology_class'] = 'Unknown'

# First subset for the top 3000 highly variable genes as determined by binomial deviance
n_top_genes = 3000
idx = adata_ref.var['binomial_deviance'].values.argsort()[-n_top_genes:]
mask = np.zeros(adata_ref.var_names.shape, dtype=bool)
mask[idx] = True
adata_ref.var['highly_deviant'] = mask
adata_ref_hvg = adata_ref[:, adata_ref.var['highly_deviant']]

# Get the intersection of query and ref var names
hvg_vars = adata_query.var_names[adata_query.var_names.isin(adata_ref_hvg.var_names)]

print(len(hvg_vars))
adata_query_hvg = adata_query[:, hvg_vars].copy()
adata_ref_hvg = adata_ref_hvg[:, hvg_vars].copy()

scvi.model.SCVI.setup_anndata(adata_ref_hvg, layer="counts")

scvi_ref = scvi.model.SCVI(
    adata_ref_hvg,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)
scvi_ref.train()

SCANVI_LABELS_KEY = "labels_scanvi"
adata_ref_hvg.obs[SCANVI_LABELS_KEY] = adata_ref_hvg.obs["cell_ontology_class"].values

scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_ref,
    unlabeled_category="Unknown",
    labels_key=SCANVI_LABELS_KEY,
)

scanvi_ref.train(max_epochs=200, n_samples_per_label=100)

scanvi_query = scvi.model.SCANVI.load_query_data(adata_query_hvg, scanvi_ref)

scanvi_query.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0},
    check_val_every_n_epoch=10,
)

SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

adata_query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()

adata_query.write(data_directory + 'reyes22_old_unspliced.h5ad', compression='gzip')

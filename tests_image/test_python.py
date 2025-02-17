#!/opt/conda/bin/python

import anndata2ri
import numpy as np
import scanpy as sc
import scipy as sp
from scipy.cluster.vq import kmeans2
from scipy.sparse import random
from sklearn.decomposition import PCA
import session_info
import spatialdata_io

import os
print('Current working directory:')
print(f'  {os.getcwd()}')

# anndata2ri
print('Test anndata2ri')
adata = sc.AnnData(X=sp.sparse.csr_matrix([[0.0, 1.0], [1.0, 0.0]]))
anndata2ri.py2rpy(adata)

# reproducibility sklearn
print('Reproducibility of pca for sparse data using sklearn')
np.random.seed(42)
X = random(2700, 32738, density=2286884/(2700*32738), format='csr', dtype='int32', data_rvs=lambda s: np.random.randint(0, 501, size=s))

pca_ = PCA(n_components=50, svd_solver='arpack', random_state=0)
X_pca_0 = pca_.fit_transform(X).copy()


for i in range(3):
    print(i)
    pca_ = PCA(n_components=50, svd_solver='arpack', random_state=0)
    X_pca = pca_.fit_transform(X).copy()
    try:
        np.testing.assert_equal(X_pca_0, X_pca)
    except AssertionError:
        session_info.show(html=False)
        raise AssertionError(f"Re-computed PCA is not the same. Maximum absolute difference: {np.abs(X_pca_0 - X_pca).max()}")


# now scanpy tests

# reproducibility dense
print('Reproducibility of pca for dense data using scanpy')

adata = sc.datasets.pbmc3k()

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)

X = adata.obsm['X_pca'].copy()


for i in range(3):
    print(i)
    sc.pp.pca(adata)
    try:
        np.testing.assert_equal(X, adata.obsm['X_pca'])
    except AssertionError:
        raise AssertionError(f"Re-computed PCA for dense is not the same. Maximum absolute difference: {np.abs(X - adata.obsm['X_pca']).max()}")

# reproducibility sparse
print('Reproducibility of pca for sparse data using scanpy')

adata = sc.datasets.pbmc3k()
#adata.X = adata.X

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)

X = adata.obsm['X_pca'].copy()

for i in range(3):
    print(i)
    sc.pp.pca(adata)
    try:
        np.testing.assert_equal(X, adata.obsm['X_pca'])
    except AssertionError:
        raise AssertionError(f"Re-computed PCA for sparse is not the same. Maximum absolute difference: {np.abs(X - adata.obsm['X_pca']).max()}")

# basic pipeline

adata = sc.datasets.blobs()

sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, 'leiden')

sc.pl.umap(adata)

# harmonypy: reproducibility
# see https://github.com/slowkow/harmonypy/issues/24
print('Reproducibility of Harmony')
adata = sc.datasets.blobs()

def cluster_fn(data, K):
    centroid, label = kmeans2(data, K, minit='++', seed=0)
    return centroid

sc.pp.pca(adata)
sc.external.pp.harmony_integrate(adata, 'blobs', cluster_fn=cluster_fn)

X_harmony = adata.obsm['X_pca_harmony'].copy()

for i in range(3):
    print(i)
    sc.external.pp.harmony_integrate(adata, 'blobs', cluster_fn=cluster_fn)
    try:
        np.testing.assert_equal(
            X_harmony, 
            adata.obsm['X_pca_harmony']
        )
    except AssertionError:
        max_diff = np.abs(X_harmony - adata.obsm['X_pca_harmony']).max()
        raise AssertionError(f"Re-computed Harmony is not the same. Maximum absolute difference: {max_diff}")
    
# Finally show the loaded packages
session_info.show(html=False)

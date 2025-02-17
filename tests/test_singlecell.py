# filepath: /Users/rost/code/singlecell-notebook/tests/test_singlecell.py
import anndata2ri
import numpy as np
import scanpy as sc
import scipy as sp
from scipy.cluster.vq import kmeans2
from scipy.sparse import random
from sklearn.decomposition import PCA
import session_info
import spatialdata_io
import pytest

def test_anndata2ri():
    """Test anndata2ri conversion."""
    adata = sc.AnnData(X=sp.sparse.csr_matrix([[0.0, 1.0], [1.0, 0.0]]))
    anndata2ri.py2rpy(adata)

def test_pca_reproducibility_sklearn():
    """Test reproducibility of PCA for sparse data using sklearn."""
    X = random(2700, 32738, density=2286884/(2700*32738), format='csr', dtype='int32', data_rvs=lambda s: np.random.randint(0, 501, size=s))
    pca_ = PCA(n_components=50, svd_solver='arpack', random_state=0)
    X_pca_0 = pca_.fit_transform(X).copy()

    for i in range(3):
        pca_ = PCA(n_components=50, svd_solver='arpack', random_state=0)
        X_pca = pca_.fit_transform(X).copy()
        np.testing.assert_equal(X_pca_0, X_pca)

def test_pca_reproducibility_scanpy_dense():
    """Test reproducibility of PCA for dense data using scanpy."""
    adata = sc.datasets.pbmc3k()
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    X = adata.obsm['X_pca'].copy()

    for i in range(3):
        sc.pp.pca(adata)
        np.testing.assert_equal(X, adata.obsm['X_pca'])

def test_pca_reproducibility_scanpy_sparse():
    """Test reproducibility of PCA for sparse data using scanpy."""
    adata = sc.datasets.pbmc3k()
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    X = adata.obsm['X_pca'].copy()

    for i in range(3):
        sc.pp.pca(adata)
        np.testing.assert_equal(X, adata.obsm['X_pca'])

def test_basic_pipeline():
    """Test basic scanpy pipeline."""
    adata = sc.datasets.blobs()
    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, 'leiden')
    sc.pl.umap(adata)

def test_harmony_reproducibility():
    """Test reproducibility of Harmony integration."""
    adata = sc.datasets.blobs()

    def cluster_fn(data, K):
        centroid, label = kmeans2(data, K, minit='++', seed=0)
        return centroid

    sc.pp.pca(adata)
    sc.external.pp.harmony_integrate(adata, 'blobs', cluster_fn=cluster_fn)
    X_harmony = adata.obsm['X_pca_harmony'].copy()

    for i in range(3):
        sc.external.pp.harmony_integrate(adata, 'blobs', cluster_fn=cluster_fn)
        np.testing.assert_equal(X_harmony, adata.obsm['X_pca_harmony'])


def test_session_info():
    """Print session information."""
    session_info.show(html=False)

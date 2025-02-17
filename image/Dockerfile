FROM quay.io/jupyter/datascience-notebook

USER root

RUN curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v1.6.40/quarto-1.6.40-linux-amd64.deb && \
    dpkg -i quarto-1.6.40-linux-amd64.deb && \
    rm quarto-1.6.40-linux-amd64.deb

USER ${NB_UID}

RUN conda config --add channels bioconda

RUN mamba install --yes \
    # Python \
    'anndata' \
    'anndata2ri' \
    'arboreto' \
    'bbknn' \
    'cellpose' \
    'celltypist' \
    'decoupler' \
    'fastcluster' \
    'gseapy' \
    'harmonypy' \ 
    'leidenalg' \
    'loompy' \
    'louvain' \
    'mscorefonts' \
    'mudata' \
    'multicore-tsne' \
    'muon' \
    'openpyxl' \
    'pybiomart' \
    'python-igraph' \
    'scanorama' \
    'scanpy' \
    'scirpy' \
    'scrublet' \
    'scvi-tools' \
    'session-info' \
    'spatialdata' \
    'spatialdata-io' \
    'spatialdata-plot' \
    'squidpy' \
    'zarr' \
    # R \
    'bioconductor-singlecellexperiment' \
    'bioconductor-monocle' \
    'bioconductor-mast' \
    'bioconductor-scater' \
    'bioconductor-scran' \
    'r-seurat' && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Update packages where their current versions have bugs
RUN pip --no-cache-dir install \
    spatial-image \
    spatialdata \
    spatialdata-io \
    spatialdata-plot \
    --upgrade
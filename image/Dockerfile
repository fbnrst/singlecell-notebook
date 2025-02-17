FROM quay.io/jupyter/datascience-notebook:latest

# USER root
# RUN curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v1.6.40/quarto-1.6.40-linux-amd64.deb && \
#     dpkg -i quarto-1.6.40-linux-amd64.deb && \
#     rm quarto-1.6.40-linux-amd64.deb

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
    'scanpy>=1.10.4' \
    'scirpy' \
    'scrublet' \
    'scvi-tools' \
    'session-info' \
    # R \
    'bioconductor-singlecellexperiment' \
    'bioconductor-monocle' \
    'bioconductor-mast' \
    'bioconductor-scater' \
    'bioconductor-scran' \
    'r-seurat>=5.1.0' && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"
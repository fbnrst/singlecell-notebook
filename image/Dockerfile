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
    'celltypist' \
    'decoupler' \
    'gseapy' \
    'harmonypy' \ 
    'leidenalg' \
    'loompy' \
    'louvain' \
    'mscorefonts' \
    'mudata' \
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

RUN pip install --no-cache-dir fastcluster

USER root

ARG CONTAINER_BUILD_DATE
ARG CONTAINER_IMAGE_NAME
ARG CONTAINER_GIT_REPOSITORY
ARG CONTAINER_GIT_COMMIT
ARG CONTAINER_VERSION

ENV CONTAINER_BUILD_DATE=${CONTAINER_BUILD_DATE}
ENV CONTAINER_IMAGE_NAME=$CONTAINER_IMAGE_NAME
ENV CONTAINER_GIT_REPOSITORY=$CONTAINER_GIT_REPOSITORY
ENV CONTAINER_GIT_COMMIT=$CONTAINER_GIT_COMMIT
ENV CONTAINER_VERSION=$CONTAINER_VERSION

LABEL org.opencontainers.image.created=$CONTAINER_BUILD_DATE \
      org.opencontainers.image.revision=$CONTAINER_GIT_COMMIT \
      org.opencontainers.image.version=$CONTAINER_VERSION

USER ${NB_UID} 
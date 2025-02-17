# singlecell-notebook

singlecell-notebook is a Jupyter Docker Stack image for single-cell sequencing analysis. It is based on the [jupyter/datascience-notebook](quay.io/jupyter/datascience-notebook) image.

To start a JupyterLab server with the latest `singlecell-notebook` container, run:

```bash
docker run -ti  --rm -p 8888:8888 -v "${PWD}":/home/jovyan/work quay.io/fbnrst/singlecell-notebook:latest
```

This will start a JupyterLab server and print a URL to access it in your browser. The current directory will be mounted as `/home/jovyan/work` in the container. See [Docker stacks documentation](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html) for more information.

<!-- For production use, rather use specific versions of the container, e.g.:

```bash
docker run -ti  --rm -p 8888:8888 -v "${PWD}":/home/jovyan/work git.mpi-cbg.de:5050/huch_lab/docker/singlecell:v0.0.1
``` -->

Find the available versions on the [releases page](https://quay.io/repository/fbnrst/singlecell-notebook).

## License
This container image is licensed under the BSD 3-Clause License.  
It based on the [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io), which are licensed under the BSD 3-Clause License.
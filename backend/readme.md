# Gauss Thinning Backend

This project contains a small websocket server which does local gauss thinning for the frontend.

## Installation

The easiest way to install the libraries is trough the [conda](https://anaconda.org/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) python package manager.

All libraries are part of the channel [conda forge](https://conda-forge.org/), which we advise to add to your conda channels by:

```bash
conda config --add channels conda-forge
```

This step allows to install any conda forge package simply with `conda install <package>`.

To install all our packages just run:

```bash
conda install igl
pip install websockets
```

## Running the backend

```bash
python src/main.py
```

## Running the backend with docker

```bash
docker build -t gauss_backend .
docker run -dp 5678:5678 gauss_backend
```

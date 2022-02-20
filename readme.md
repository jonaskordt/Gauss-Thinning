# Algorithmic folding course - In-depth assignment on Gauss Thinning

## Project

I plan to implement local Gauss-Thinning in Python together with a frontend for user input.

## Structure

The repo contains some jupyter notebooks which can execute global and local Gauss Thinning on example meshes.
Additionally, the folder `frontend` contains code for a React frontend which uses Three.js to render the meshes and the folder `backend` contains a small python websocket server which will do the local Gauss Thinning triggered by the frontend in the future.

## Resources

- [Gauss Thinning project website](https://igl.ethz.ch/projects/gauss-thinning/)
- [Gauss Thinning paper](https://igl.ethz.ch/projects/gauss-thinning/GaussThinning_Paper.pdf)
- [Gauss Thinning paper presentation](https://www.youtube.com/watch?v=k0RVs_FKYd4)
- [Gauss Thinning repository (C++)](https://github.com/FloorVerhoeven/DevelopableApproximationViaGaussImageThinning)

## Installation (notebooks)

The easiest way to install the libraries is trough the [conda](https://anaconda.org/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) python package manager.

All libraries are part of the channel [conda forge](https://conda-forge.org/), which we advise to add to your conda channels by:

```bash
conda config --add channels conda-forge
```

This step allows to install any conda forge package simply with `conda install <package>`.

To install all our packages just run:

```bash
conda install igl
conda install meshplot
conda install jupyter
```

If you have problems with the pythreejs rendering, check out [their documentation](https://github.com/jupyter-widgets/pythreejs).

Only Chrome is supported right now.
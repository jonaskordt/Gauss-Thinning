# Algorithmic folding course - In-depth assignment on Gauss Thinning
=======

## Project

I plan to implement Gauss-Thinning in Python and then see where to go from there.

## Resources

- [Gauss Thinning project website](https://igl.ethz.ch/projects/gauss-thinning/)
- [Gauss Thinning paper](https://igl.ethz.ch/projects/gauss-thinning/GaussThinning_Paper.pdf)
- [Gauss Thinning paper presentation](https://www.youtube.com/watch?v=k0RVs_FKYd4)
- [Gauss Thinning repository (C++)](https://github.com/FloorVerhoeven/DevelopableApproximationViaGaussImageThinning)

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
conda install meshplot
conda install jupyter
```

If you have problems with the pythreejs rendering, check out [their documentation](https://github.com/jupyter-widgets/pythreejs).

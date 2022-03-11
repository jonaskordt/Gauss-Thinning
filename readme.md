# Algorithmic folding course - In-depth assignment on Gauss Thinning

## Project

We implemented Local Gauss Thinning and an interactive interface.

## Structure

The repo contains some jupyter notebooks which can execute global and local Gauss Thinning on example meshes.
Additionally, the folder `frontend` contains code for a React frontend which uses Three.js to render the meshes and the folder `backend` contains a small python websocket server which will do the local Gauss Thinning triggered by the frontend in the future.

## Resources

- [Gauss Thinning project website](https://igl.ethz.ch/projects/gauss-thinning/)
- [Gauss Thinning paper](https://igl.ethz.ch/projects/gauss-thinning/GaussThinning_Paper.pdf)
- [Gauss Thinning paper presentation](https://www.youtube.com/watch?v=k0RVs_FKYd4)
- [Gauss Thinning repository (C++)](https://github.com/FloorVerhoeven/DevelopableApproximationViaGaussImageThinning)

## Compatibility

Only Chrome and Safari are supported right now.

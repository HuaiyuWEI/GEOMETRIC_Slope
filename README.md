# Overview
This repository contains the modified MITgcm source code used in the manuscript "Parameterizing eddy buoyancy fluxes across prograde shelf/slope fronts using a slope-aware GEOMETRIC closure". The code is based on the GEOMETRIC parameterization ([David et al., 2012][https://www.sciencedirect.com/science/article/pii/S1463500319301775?via%3Dihub)]) implemented into MITgcm by [Julian Mak](https://github.com/julianmak/GEOMETRIC_code). This code adds a "slope-aware" modification to the GEOMETRIC parameterization to improve the prediction of eddy buoyancy diffusivity over steep continental slopes ([Wei et al., 2022](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022MS003229)).  

# Dependency
MITgcm (version [checkpoint68n](https://zenodo.org/record/762177))

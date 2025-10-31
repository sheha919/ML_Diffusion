# Predicting OH<sup>-</sup> Diffusion in the dual-cation AEM using ML 
Author: Shehani Wetthasinghe

Last modified: 10/31/2025

![cations](https://github.com/sheha919/ML_Diffusion/blob/main/img/dual_cat_sys.png)

This repository contains Python scripts to generate the systems/input files for Molecular Dynamics simulation in DFTB+ and conduct the analysis/further calculations based on the performed simulations.

* [gen_graphenebox/generate_sheets.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/generate_sheets.ipynb): Create the initial system (two graphene sheets without cations) 
* [gen_graphenebox/add_cations.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/add_cations.ipynb)    : Add cations to the system
* [gen_files/gen_mdsim.py](https://github.com/sheha919/ML_Diffusion/blob/main/gen_files/gen_mdsim.py)               : Create initial and extended input files
* [gen_data/chem_prop.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_data/chem_prop.ipynb)             : Generate chemical descriptors using RDkit 

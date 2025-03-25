# Predicting OH<sup>-</sup> Diffusion in the dual-cation AEM using ML 
Author: Shehani Wetthasinghe

Last modified: 03/17/2025

![cations](https://github.com/sheha919/ML_Diffusion/blob/main/gen_data/mol_str/cat_grid.png)

Here we are interested on 40 cations that can be used in Anion Exchange Membranes(AEM) of Alkaline Fuel Cell (AFC) to increase the diffusion of OH<sup>-</sup> within the membrane.

I have used following scripts to create the systems/input files for Molecular Dynamics simulation in DFTB+ and conduct the analysis/further calculations based on the performed simulations.

* [gen_graphenebox/generate_sheets.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/generate_sheets.ipynb): Create the initial system (two graphene sheets without cations) 
* [gen_graphenebox/add_cations.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/add_cations.ipynb)    : Add cations to the system
* [gen_files/gen_mdsim.py](https://github.com/sheha919/ML_Diffusion/blob/main/gen_files/gen_mdsim.py)               : Create initial and extended input files
* [gen_data/chem_prop.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_data/chem_prop.ipynb)             : Generate chemical descriptors using RDkit 

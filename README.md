# Predicting OH<sup>-</sup> Diffusion in the dual-cation AEM using ML 
Author: Shehani Wetthasinghe

Last modified: 10/31/2025

![cations](https://github.com/sheha919/ML_Diffusion/blob/main/img/dual_cat_sys.png)

This repository contains Python scripts to generate the systems/input files for Molecular Dynamics simulation in DFTB+ and conduct the analysis/further calculations based on the performed simulations.

Steps to generate the systems/input files:
1. Generate the initial configuration containing two graphene sheets with a specified separation along the Z-axis and defined simulation box dimensions in the X and Y directions. This script also calculates the number of water molecules required to achieve the specified water number density for a given box size with the selected single or dual cations: [gen_graphenebox/generate_sheets.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/generate_sheets.ipynb)
2. Translate the system to the center using IQmol
3. Add selected single or dual cations with the preferred orientation.
     * single cation to the top graphene sheet: [gen_graphenebox/add_cations.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/add_cations.ipynb)
     * dual cations to the top graphene sheet: [gen_graphenebox/add_dual_cations.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/add_dual_cations.ipynb)
     * one cation to the top graphene sheet and the second one to the bottom graphene sheet:[gen_graphenebox/add_dual_cations_inv.ipynb](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/add_dual_cations_inv.ipynb)
4. Add the calculated number of water molecules (given in step1) to the system using Spartan
5. Generate an initial DFTB+ input file containing frozen graphene sheets and cation(s) to facilitate the construction of multiple systems with distinct water molecule orientations: gen_infile_frozen.py[gen_graphenebox/gen_infile_frozen.py](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/gen_infile_frozen.py)
6. Generate multiple input files for thermalization simultaneously using the distinct water molecule orientations from step 5: [gen_graphenebox/gen_sim_therm.py](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/gen_sim_therm.py)
7. Analyze the thermalization statistics to confirm that the system has reached proper thermal equilibrium: [gen_graphenebox/therm_analysis.py](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/therm_analysis.py)
8. Generate multiple input files for molecular dynamics simulations simultaneously using the final geometry extracted from the thermalization step: [gen_graphenebox/gen_sim.py](https://github.com/sheha919/ML_Diffusion/blob/main/gen_graphenebox/gen_sim.py)
9. Check for thermalization statistics at the end of dynamics simulations.

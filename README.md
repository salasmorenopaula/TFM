# Heavy Dark Matter Analysis with Gammapy
Author: Paula Salas Moreno (pausalas@ucm.es)
This repository contains the code developed for the Master's Thesis focused on indirect detection of heavy dark matter in the Draco dwarf galaxy using the gammapy package and an Effective Field Theory framework.

## Contents

The repository includes the following key files:

### 1. Heavy_DM_model.ipynb
A Jupyter Notebook containing the majority of the code used in the thesis. It is divided into two main parts:

- **Part 1:** 
  - Calculation of the J-factor for Draco, Ursa Mayor II and Wilman I.
  - Point source analysis.
  - Draco observation simulation with CTA for a 10 TeV dark matter particle.
  - Generation of significance and counts maps.

- **Part 2:**  
  - Analysis of annihilation into gamma-ray lines.
  - Computation of parameters relevant to the freeze-out epoch.
  - Calculation and plotting of the annihilation cross section ⟨σv⟩ as a function of WIMP mass.
  - Comparison with experimental limits from Fermi-LAT, H.E.S.S., MAGIC, and CTA.
  - Computation of the dilusion factor.

### 2. simulation_draco_multiple_masses.py
A Python script that performs 30 simulations of gamma-ray observations of Draco for a user-defined dark matter mass (default channel: annihilation into two Higgs bosons, but this can be changed).  
It computes the significance and the upper limit on ⟨σv⟩ for each simulation and saves the results to a CSV file.

### 3. Exclusion_curve_sigma_Draco_simulation_for_several_masses.ipynb
A Jupyter Notebook that uses the results from the multiple DM mass simulations to:
- Plot the exclusion curve (⟨σv⟩ vs DM mass).
- Create a summary table with significance values and upper limits for each mass.
  
## Installation and Requirements

To run the code, a standard scientific Python environment is needed. The following packages are required:

- `numpy`
- `matplotlib`
- `pandas`
- `scipy`
- `astropy`
- `gammapy`
- `math` (built-in)
- `sympy`
- `argparse` (built-in, used in the `.py` script)
  
##  Data Requirements

### 1. Gammapy Data

Some simulations rely on template data provided by the Gammapy team.

Download the dataset from the official GitHub repository:  
🔗 [https://github.com/gammapy/gammapy-data](https://github.com/gammapy/gammapy-data)

Although only the `dark_matter_spectra` folder may be strictly necessary, downloading the full repository is recommended.

After downloading, set the GAMMAPY_DATA environment variable to the folder path.  

### 2. CTA Instrument Response Functions (IRFs)

To perform CTA simulations, the IRFs must be downloaded from the official Zenodo record:  
🔗 [https://zenodo.org/records/5499840](https://zenodo.org/records/5499840)

Make sure that your notebooks and scripts correctly reference the location of the downloaded IRF files.

---

 The code is designed to be modular and reproducible.  
All results are generated using gammapy and are intended for use in indirect dark matter searches via gamma-ray observations.


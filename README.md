# LSS-bispectrum

## Bispectrum Calculation and Analysis in Redshift Space  

This repository contains scripts and analysis related to the calculation of the bispectrum in redshift space (RSD) for large-scale structure (LSS) studies. The work is based on **Quijote simulations** and investigates the bispectrum across various triangle configurations, focusing on the impact of redshift-space distortions.  

## Contents  
- **`bispec_50QSims.py`**: Python script to compute the bispectrum in redshift or real space.  
- **`Bk configuration analysis/ `**: Jupyter notebooks analyzing bispectrum results for different triangle configurations along with plots
- **`data/`**: Sample simulation data from Quijote Simulations.  


## Features  
- **Redshift-Space Distortion (RSD) Modeling**: Leverages tools to apply RSD corrections and compute the bispectrum.  
- **Triangle Configuration Analysis**: Examines bispectrum differences across equilateral, squeezed, and flat configurations.  


## Prerequisites  
- **Python 3.8+**  
- Required packages:  
  - `numpy`  
  - `pandas`  
  - `matplotlib`  
  - `Pylians3`  

> **Note**: Pylians3 is designed to work on macOS and Linux. If you're using Windows, it's recommended to run the scripts via **Windows Subsystem for Linux (WSL)** for compatibility.



## Results  
Key findings from the analysis include:  
- **Sensitivity to RSD**: The bispectrum changes significantly in redshift space compared to real space.  
- **Configuration Dependence**: Squeezed triangle configurations show greater sensitivity to non-linear effects.  

## Acknowledgments  
This work was conducted as part of my Master's dissertation and I would like to give thanks to Chris Clarkson and Chris Addis for their guidance. The Quijote simulations were used as the primary data source. Special thanks to the developers of Pylians3 for providing essential tools for large-scale structure analysis


# Optimizing Dementia Risk Reduction

This repository contains the analysis code for the paper Optimizing Dementia Risk Reduction: Balancing Sleep and Physical Activity Trade-offs by Yiallourou et al. This project made use of the UK Biobank resource. Data analysis was performed by Beaudan Campbell-Brown and Lachlan Cribb.

## Description

Code for executing the analysis and providing key outputs are contained in R files with the _bootstrap suffix (e.g., dem_bootstrap.R for the primary dementia risk analysis). Helper functions for these files will be contained in the respective files with _models suffix and in the utils.R script. Create_dataset.R and prepare_data.R are files for selecting and preparing relevant variables from the UK Biobank resource for use with this project.

## Packages

All packages and package dependencies used in the project are available in the renv.lock. See the R [renv](https://rstudio.github.io/renv/articles/renv.html) package for details.
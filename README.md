# Comparison of methods for adjusting for exposure-affected time-varying confounding under misspecification – Monte-Carlo Simulations
## Stata and R Analysis Code

This repository contains the Stata and R code used in the misspecification simulation by Clare et al. 2019

The Stata code creates a series of quasi-random datasets using a pre-specified data structure.
Analysis code runs all analyses on those datasets, and saves the results. Note that the code is written to run on Google Compute clusters, using a Linux OS (in order to run the syntax on a Windows-based machine, some changes to the way parallel processing is required (because Windows is not compatible with 'FORK').

| Description | Code |
| --- | --- |
| S1 - Data creation - Stata Code | [Data creation code](Code/S1_data_creation.do) |
| S2 - Analysis - R Code | [Analysis code](Code/S2_analysis_code.R) |




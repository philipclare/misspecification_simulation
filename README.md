# Comparison of methods of adjusting for time-varying confounding under misspecification – A Monte-Carlo simulation study
## Stata and R Analysis Code

This repository contains the Stata and R code used in the misspecification simulation by Clare et al. 2019

The Stata code creates a series of quasi-random datasets using a pre-specified data structure.
Analysis code runs all analyses on those datasets, and saves the results. Note that the code is written to run on Google Compute clusters, using a Linux OS (in order to run the syntax on a Windows-based machine, some changes to the way parallel processing is required (because Windows is not compatible with 'FORK').

Two types of standard error estimates were used, so two sets of analysis code are included. The first calculates standard errors using bootstrapping. The second calculates model-based standard errors, using influence curves for TMLE.

| Description | Code |
| --- | --- |
| S1 - Data creation Stata Code | [Data creation code](Code/S1_data_creation.do) |
| S2 - Analysis with bootstrap SEs - R Code | [Analysis code - Bootstrap](Code/S2_analysis_code_bootstrap.R) |
| S3 - Analysis with model-based/influence curve SEs - R Code | [Analysis code - Alternative](Code/S3_analysis_code_IC.R) |



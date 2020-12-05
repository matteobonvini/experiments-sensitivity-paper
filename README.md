# Experiments-for-Sensitivity-Paper
Simulations and data analysis for the paper "Sensitivity Analysis via the Proportion of Unmeasured Confounding." See the manuscript and the supplementary materials for further details and references. 

To run these scripts, one needs to install the "sensitivitypuc" package from Github repository "matteobonvini/sensitivitypuc" as

# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # if warnings are causing erros and they are not important one may use this
devtools::install_github("matteobonvini/sensitivitypuc")

- Section 4.1 (Simulation study):
	- 1. Run "simulation_get_truth.R". It creates the true values used in the simulation;
	- 2. Run "simulation.R". It runs the actual simulation. 

- Section 4.2 and Section 6 in Supplementary materials (Data analysis):
	- 1. Run "data analysis.R". It generates the bounds and confindence bands using our model.
	- 2. Run "data_analysis_plot.R". It generates the plots contained in the paper.
	- 3. Run "data_analysis_comparison_cinelli_hazlett2020.R". It generates the plots used in Section 6.1 from the sensitivity model in Cinelli & Hazlett (2020). 

- Section 7 in Supplementary materials (Simulations regarding power):
	- 1. Run "design_sensitivity.R". It generates the plots contained in the section. 

The folder also contains the file "plot_theme.R" that specifies the main ggplot2 theme used in all the plots; it will be sourced in all the R scripts that generate plots. The folder also contains the file "example_partition.R" that is used to create Figure 1 in the main text: it creates the scatter plots to which we have manually added the red regions. This last file can be ignored.

################################################################################################################
#                                                                                                               
# Files to reproduce the analyses of the article:                                                                  
#                                                                                                               
# Diversity mediates the responses of invertebrate density to annual drying duration and frequency              
#                                                                                                               
# Rebeca Arias-Real, Cayetano Gutiérez-Cánovas, Margarita Menéndez, Verónica Granados 
# & Isabel MuÃ±oz
#                                                                                                               
# Code written by Rebeca Arias-Real and Cayetano Gutiérrez-Cánovas                                              
# email for queries: rebeca.arias.real@gmail.com or cayeguti@um.es                                              
################################################################################################################

0_FD_functions.R
Functions to estimate Functional Diversity metrics

0_quality_funct_space_fromdist.R
R function for computing the quality of functional dendrogram and multidimensional functional spaces

1_metric_calculations.R
R script to calculate taxonomic and trait-based metrics for invertebrate communities

2_lineal_models.R
R script to run models evaluating the effect and importance of drying aspects on inverebrate communities 

3_SEM.R
R script to run SEM to explore direct and indirect impacts of drying aspects on invertebrate density

4_model_plots.R
R script to plot invertebrate metric responses to drying aspects

env.txt
Environmental variables, including hydrological variables

env_description.txt
Full description of environmental variables

final_set.txt
Dataset used in linear models and SEMs

inv_kick.txt
Invertebrate taxa from kick samples

inv_surber.txt
Invertebrate taxa from surber samples

traits.txt
Functional trait data

traits_description.txt
Description of functional traits

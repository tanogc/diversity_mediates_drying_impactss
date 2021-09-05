

# Code to reproduce the statistical analyses of "Diversity mediates the responses of invertebrate density to annual drying duration and frequency"

In this paper, using aquatic invertebrates from 33 streams across a flow-intermittence gradient, we assessed how annual (drying duration and frequency) and recent drying characteristics (duration of the last dry period and flowing duration since the last rewetting) affect the density and diversity metrics of communities and trophic groups while controlling for other key abiotic factors (dissolved oxygen and altitude). We characterized invertebrate communities using taxonomy and functional traits to capture biological features that increase vulnerability to drying. In addition, using structural equation modelling (SEM), we evaluated pathways by which drying characteristics directly impact invertebrate density and whether diversity indirectly mediates such relationships.   

This code re-creates analysis of the field dataset

## Original article

Please, use this citation to reference the database:
```
Arias-Real, R., Gutiérrez-Cánovas, C., , Menéndez, M., Granados, V. & Muñoz, I.
Diversity mediates the responses of invertebrate density to annual drying duration and frequency. 
Oikos (accepted)
```

# R files description

*	0_FD_functions.R: Functions to estimate Functional Diversity metrics
* 0_quality_funct_space_fromdist.R: R function for computing the quality of functional dendrogram and multidimensional functional spaces
* 1_metric_calculations.R: R script to calculate taxonomic and trait-based metrics for invertebrate communities
* 2_lineal_models.R: R script to run models evaluating the effect and importance of drying aspects on inverebrate communities 
* 3_SEM.R: R script to run SEM to explore direct and indirect impacts of drying aspects on invertebrate density
* 4_model_plots.R: R script to plot invertebrate metric responses to drying aspects

# Data
* env.txt: Environmental variables, including hydrological variables
* env_description.txt: Full description of environmental variables
* final_set.txt: Dataset used in linear models and SEMs
* inv_kick.txt: Invertebrate taxa from kick samples
* inv_surber.txt: Invertebrate taxa from surber samples
* traits.txt: Functional trait data
* traits_description.txt: Description of functional traits

```
Please, send questions or problems related with the use of this code to Rebeca Arias-Real (rebeca.arias.real@gmail.com) 
or Cayetano Gutiérrez-Cánovas (cayeguti@um.es).


**Code repository for "Enriched Genetic Pathways Shape Generative Structural Brain Mechanisms and Cognitive Ability".**

Alicja Monaghan [1,2], Danyal Akarca [1] & Duncan Astle [1,3]
1. MRC Cognition and Brain Sciences Unit, University of Cambridge, Cambridge, UK
2. The Alan Turing Institute, London, UK
3. Department of Psychiatry, University of Cambridge, Cambridge, UK

For any questions about the use of this repository, please contact Alicja.Monaghan@mrc-cbu.cam.ac.uk.

**Requirements**

The following installations are required to use all the attached scripts. 
* MATLAB 2022a (installation: https://uk.mathworks.com/help/install/install-products.html)
* RStudio R 2022.12.0.353 (installation: https://rstudio.com/products/rstudio/download/)
* Brain Connectivity Toolbox 2019 (installation: https://sites.google.com/site/bctnet/)
* Cytoscape v3.9.1 (installation: https://cytoscape.org/download.html)

**Data Availability Statement**

This work uses neuroimaging, genetic, and cognitive data from 2154 children from the baseline time point of the Adolescent Brain Cognitive Development 
(ABCD) study, Release 4 (http://dx.doi.org/10.15154/1523041). Access to ABCD can be requested here: https://nda.nih.gov/nda/access-data-info.html. Since ABCD is 
an open-access initiative, to reproduce our results, we provide the following derived data: structural connectomes generated by deterministic tractography 
in QSIprep 0.15.3 (Cieslak et al., 2021), framewise displacement from QSIprep, age in months, sex, site ID (recoded into random letters), polygenic scores for 
cognitive ability (calculated using a tutorial from Choi and colleagues, 2020, Nature Protocols, 15, 2759-2772), general g factor of cognitive 
ability (derived from a principal component analysis of 7 NIH Toolbox cognitive assessements of working memory, attention, executive function, and language). 
Note that all data from ABCD is anonymised, and no subject IDs are included. We've also undertaken additional anonymisation steps, such as recoding 
scanner site. 

**/Scripts**
1. threshold_and_binarise.m
2. run_and_evaluate_group_gnm.m
3. run_individual_neighbours_structural_gnm_grid_search.m
4. examine_individual_neighbours_structural_gnm_outputs.m
5. run_pls_regression_ahba_nodal_parameters.m
6. sorting_pls_regression_ahba_nodal_parameters.m
7. linking_ahba_genetics_with_gnm_parameters_and_cognition.m
8. abcd_gene_ontology_analysis.R
9. gene_ontology_and_comparative_pathway_enrichment.cys --> Note, this contains two annotation sets. The original contains the original node names/groupings, whilst the 
second contains a curated version in the manuscript. 
10. abcd_group_and_individual_gnm_visualisation.R
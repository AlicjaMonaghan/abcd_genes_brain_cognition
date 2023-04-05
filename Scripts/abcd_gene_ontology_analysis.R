# This script details gene ontology (GO) analyses for Allen Human Brain Atlas 
# (AHBA) genes predictive of mean parameterised nodal wiring costs and values, 
# respectively, alongside single-nucleotide polymorphisms included in the 
# polygenic scores (PGSs) for cognitive ability, across 1461 children from the 
# baseline Adolescent Brain Cognitive Development (ABCD) study. Note that we 
# restricted our analyses to the left hemisphere, for which all 6 AHBA donors
# provided data for, and therefore had better spatial coverage. Correspondence
# to Alicja Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk.

# STEPS: 
# 1. Setting up the work space.
# 2. Loading data.
# 3. GO analyses for parameterised nodal wiring costs.
# 4. GO analyses for parameterised nodal wiring value.
# 5. GO analyses for SNPs included in PGS for cognitive ability.
# 6. Multi-query comparative pathway enrichment for lists in steps 3-5. Note 
# that step 6 is visualised in the accompanting Cytoscape session.

### STEP 1 - Setting Up the Work Space ####
rm(list = ls())
working_directory = "abcd_genomic_variation_structural_generative_mechanisms_open/"
# SET TO WHERE YOU SAVED THIS DIRECTORY!
setwd(working_directory)

# Load the packages!
library(gprofiler2)
library(tibble)
library(purrr)

### STEP 2 - Load Data ####
# Load the summary table for the partial least square (PLS) regressions of AHBA
# genes from the left hemisphere as predictors of parameterised nodal wiring 
# costs (averaged across participants).
nodal_wiring_costs_genes = read.csv('data/nodal_wiring_costs_genes.csv',header = T)
# Load the equivalent, but for parameterized nodal wiring value.
nodal_wiring_value_genes = read.csv('data/nodal_wiring_value_genes.csv',header=T)
# Load the SNPs included in the PGS for cognitive ability, ranked by predictive
# power i.e. decreasing beta.
pgs_ranked_snps = read.table('data/pgs_ranked_snps.txt',header = T)

### STEP 3 - GO Analyses for Parameterised Nodal Wiring Costs ####
nodal_wiring_costs_gost = gost(nodal_wiring_costs_genes$ahba_gene_name, 
                               organism = "hsapiens", 
                               ordered_query = T, significant = T, 
                               user_threshold = .05, sources = "GO",
                               exclude_iea = T)
# Format the results appropriately for import into Cytoscape
nodal_wiring_costs_gost_formatted_results = 
  tibble(GO.ID = nodal_wiring_costs_gost$result$term_id,
         Description = nodal_wiring_costs_gost$result$term_name,
         p.val = nodal_wiring_costs_gost$result$p_value,
         Phenotype = "+1",
         Source = factor(nodal_wiring_costs_gost$result$source))
### STEP 4 - GO Analyses for Parameterised Nodal Wiring Value ####
nodal_wiring_value_gost = gost(nodal_wiring_value_genes$ahba_gene_name, 
                               organism = "hsapiens", 
                               ordered_query = T, significant = T, 
                               user_threshold = .05, sources = "GO",
                               exclude_iea = T)
# Format the results appropriately for import into Cytoscape
nodal_wiring_value_gost_formatted_results = 
  tibble(GO.ID = nodal_wiring_value_gost$result$term_id,
         Description = nodal_wiring_value_gost$result$term_name,
         p.val = nodal_wiring_value_gost$result$p_value,
         Phenotype = "+1",
         Source = factor(nodal_wiring_value_gost$result$source))




### STEP 5 - GO Analyses for SNPs included in PGS for Cognitive Ability ####
# Note, this GO analysis takes particularly long to run (up to 10 minutes)
abcd_cognitive_ability_gost = 
  gost(pgs_ranked_snps$SNP, ordered_query = T,
       significant = T, user_threshold = .05,
       sources = "GO", organism = "hsapiens", exclude_iea = T)

# Format the results appropriately for import into Cytoscape
abcd_cognitive_ability_gost_formatted_result = 
  tibble(GO.ID = abcd_cognitive_ability_gost$result$term_id,
         Description = abcd_cognitive_ability_gost$result$term_name,
         p.val = abcd_cognitive_ability_gost$result$p_value,
         Phenotype = "+1",
         Source = factor(abcd_cognitive_ability_gost$result$source))

### STEP 6 - Multi-Query Comparative Pathway Enrichment ####
# We will now perform a comparative enrichment of the 3 gene lists by submitting
# the gene lists to gProfiler. Note that when we must separate each list by a >QUERY_TITLE.
cost_value_cognition_pgs_comparison_list = 
  list("nodal_wiring_costs_ordered_genes" = c(nodal_wiring_costs_genes$ahba_gene_name),
       "nodal_wiring_value_ordered_genes" = c(nodal_wiring_value_genes$ahba_gene_name),
       "abcd_pgs_snps_ordered" = c(pgs_ranked_snps$SNP))

cost_value_cognition_pgs_comparison_gost = 
  gost(query = cost_value_cognition_pgs_comparison_list,
       organism = "hsapiens", ordered_query = T, multi_query = T,
       significant = T, exclude_iea = T, user_threshold = .01,
       sources = "GO")

# Extract the p-values. This produces 3 sets of p-values, corresponding to 
# testing for significant differences in enrichment between the gene lists.
cost_value_cognition_pgs_comparison_gost_pval = 
  transpose(cost_value_cognition_pgs_comparison_gost$result$p_values, 
            c("Cost", "Value", "Cognition")) %>% map_dfr(unlist)
# Format nicely. This is Supplementary Table 9.
cost_value_cognition_pgs_comparison_gost_formatted = 
  tibble(GO.ID = cost_value_cognition_pgs_comparison_gost$result$term_id,
         Description = cost_value_cognition_pgs_comparison_gost$result$term_name,
         p.val.Cost = c(cost_value_cognition_pgs_comparison_gost_pval[,1])$Cost,
         p.val.Value = c(cost_value_cognition_pgs_comparison_gost_pval[,2])$Value,
         p.val.Cognition = c(cost_value_cognition_pgs_comparison_gost_pval[,3])$Cognition,
         Phenotype = "+1",
         Source = factor(cost_value_cognition_pgs_comparison_gost$result$source))
# The comparative enrichment is included as a separate Cytoscape session. It 
# takes 3 inputs - the p-values associated with parameterised nodal wiring costs,
# value, and the SNPs included in the PGS for cognitive ability. 

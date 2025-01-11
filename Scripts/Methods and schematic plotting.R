library(ggplot2)
library(paletteer)
library(ggpointdensity)
library(R.matlab)
library(igraph)
library(ggseg)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggsegSchaefer)
library(ggtext)
library(ggrepel)
library(ggpubr)
library(tidyverse)
rm(list=ls())

setwd('/Users/alicjamonaghan/Desktop/abcd_genes_brain_cognition/')
# Load the polygenic scores for cognitive ability across the whole ABCD sample
abcd_pgs_df = read.table('data/ABCD_cognitive_ability_pgs.txt', header = TRUE)
# Load the generative network parameters and polygenic scores
gnm_pgs_df = read.csv('data/gnm_parameters_pgs_covariates.csv')
# Load the seed network used for all simulations, and extract the Schaefer 100-
# node variant. 
seed = readMat('data/seed_across_parcellations.mat')[["seed"]][[1]]
# And load the consensus network (group simulation target)
consensus = readMat('data/group_target_across_parcellations.mat')[["target"]][[1]]
# Find the coordinates for the Schaefer 100-node parcellation
schaefer100_coords = readMat('data/schaefer100x17_1mm_info.mat')[["schaefer100x17.1mm.info"]]
coords = data.frame(x = schaefer100_coords[[2]], y = schaefer100_coords[[3]], z = schaefer100_coords[[4]])
# Load the optimal participant-level eta and gamma estimates
eta_gamma_participant_level = read.csv('data/optimal_participant_level_gnm_parameters.txt') %>%
  # Convert into a long format and rename columns
  melt() %>% rename("parameter" = variable) %>% mutate(parameter = factor(parameter))
# Load the nodal AHBA expression data
ahba_expression = read.csv('data/AHBA_expression_rnaseq_schaefer100_cleaned.csv')
# Load the functional enrichment results for parameterized nodal wiring costs
nodal_costs_func_enrichment = read.csv(
  'data/nodal_wiring_costs_gost_formatted_results.csv')
# And the equivalent for parameterized nodal wiring value
nodal_value_func_enrichment = read.csv(
  'data/nodal_wiring_value_gost_formatted_results.csv'
)

### Part 1 - Polygenic scores for general cognitive ability ####
# Multiply the polygenic scores by 10000 for plotting the density
abcd_pgs_df$SCORE = abcd_pgs_df$SCORE*10000
pgs_distribution_plot = ggplot(
  data=abcd_pgs_df, mapping=aes(x = SCORE, y = PHENO)) + geom_pointdensity() + 
  labs(x = bquote("Polygenic scores for cognitive ability "(x10^-4)),
       y = "Cognitive ability") +
  theme(panel.background = element_blank(), legend.position = "none",
        axis.text = element_text(size=15), axis.title = element_text(size = 15),
        axis.line = element_line(color="black")) +
  scale_color_paletteer_c("pals::ocean.matter")
ggsave(filename = "pgs_distribution_plot.png", plot = pgs_distribution_plot,
       height = 5, width = 5, dpi = 700)

### Part 2 - Distribution of participant-level eta and gamma parameters ####
# Create two box plots to show the distribution of eta and gamma...
ggplot(eta_gamma_participant_level, mapping = aes(x = value)) +
  facet_grid(~parameter, scales="free") + geom_boxplot()


### Part 3 - Seed and consensus networks ####
seed_graph = graph_from_adjacency_matrix(seed+consensus, 'undirected', weighted = TRUE)
# Highlight the seed network amidst the consensus network
node_values = rowSums(seed+consensus)
node_colors = ifelse(node_values >= 2, "blue", "pink")
V(seed_graph)$color = node_colors
plot(seed_graph, vertex.color = V(seed_graph)$color, vertex.label = NA, layout = data.matrix(coords),
     vertex.size = 4)

V(seed_graph)$color = ifelse(V(seed_graph)$weight > 1, 'grey', 'pink')

g2 = simplify(seed_graph)
E(seed_graph)$color[E(seed_graph)$weight > 1] = 'darkgreen'
E(seed_graph)$color[E(seed_graph)$weight == 1] = 'grey'
plot(seed_graph, vertex.label = NA, axes = FALSE, vertex.size = 3)



plot(seed_graph, vertex.label.color="black", vertex.label = NA,
     
     vertex.color=c( "pink", "skyblue")[1+(V(seed_graph)$weight == 1)] ) 




seed_plot = plot(
  seed_graph, vertex.size = 0, vertex.shape = 'none', 
  layout=data.matrix(coords), margin=-1, vertex.label=NA, axes = FALSE)

### Part 4 - Visualize nodal distribution of AHBA regions previously linked to structural brain development ####
# Add the labels for the Schaefer 100-node parcellation
ahba_expression$region = unlist(schaefer100_coords[[1]])
ahba_genes_to_visualise = c("SEMA3A", "CCDC88C", "NUAK1")
for (gene_idx in 1:length(ahba_genes_to_visualise)){
  gene_to_visualise = ahba_genes_to_visualise[gene_idx]
  # Extract the expression data for the gene of interest and the region label
  gene_idx_in_ahba = which(colnames(ahba_expression) %in% gene_to_visualise)
  region_label_idx = which(colnames(ahba_expression) %in% "region")
  expression_df = ahba_expression[,c(gene_idx_in_ahba, region_label_idx)]
  # Visualize the gene's distribution in the left hemisphere
  ahba_nodal_gene_expression_visualisation = ggseg(
    .data = expression_df, atlas = schaefer17_100, 
    mapping = aes(fill = get(gene_to_visualise)), hemisphere = "left") + 
    theme_void() + labs(fill="") + 
    paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1) +
    theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))
  ggsave(filename = sprintf("%s.nodal.gene.expression.png", gene_to_visualise),
         height = 2.5, width = 2.5, dpi = 700, plot = ahba_nodal_gene_expression_visualisation)
}

### Part 5 - Visualize gene enrichment of nodal wiring parameters ####
# Plot each term as an individual circle whose diameter is proportional to the
# term size. The inside will be colored according to gene ontology source.
parameterised_wiring_terms = c("costs", "value")
# Create a list to hold the plots
parameterised_wiring_terms_plot_list = vector("list", 2)
for (i in 1:length(parameterised_wiring_terms)){
  # Retrieve the data frame to plot
  wiring_df = get(paste0("nodal_", parameterised_wiring_terms[i], "_func_enrichment")) %>%
    # Rename factor levels for plotting
    mutate('Gene ontology category' = fct_recode(
      Source, "Biological processes" = "GO:BP",  "Cellular components" = "GO:CC",
      "Molecular functions" = "GO:MF"))
  # Customize the title according to which parameter is plotted and the number
  # of functional enrichment categories
  if (parameterised_wiring_terms[i] == "costs"){
    plot_title = bquote(.(nrow(wiring_df))~"GO functional enrichments for"~eta)
  } else{
    plot_title = bquote(.(nrow(wiring_df))~"GO functional enrichments for"~gamma)
  }
  # Plot the functional enrichment...
  parameterised_wiring_terms_plot_list[[i]] = 
    ggplot(data=wiring_df, aes(x = p.val, y = GeneRatio)) +
    geom_point(aes(size=TermSize, colour = Source, fill = Source)) +
    labs(x = "Enrichment *p*-value", y = "Gene ratio",
         title = plot_title) +
    # "Least significant" results will be on the left-hand side
    scale_x_reverse() + 
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"),
          axis.title = ggtext::element_markdown(size=15), axis.text = element_text(size=15),
          plot.title = element_text(size=15, hjust=0.5)) +
    scale_colour_manual(values = c("GO:BP" = "#2F0F3E", "GO:CC" = "#CE4356", "GO:MF" = "#F29567")) +
    # Label the most important hits i.e. those with the smallest 
    geom_text_repel(data = subset(wiring_df, 0.25 < GeneRatio & GeneRatio < 1), aes(label=Description),
                    box.padding = 0.5, point.padding = 0.5, min.segment.length = 0)
}
# Combine the two plots together!
parameterised_gene_ontology_func_enrichment_plot = ggarrange(
  parameterised_wiring_terms_plot_list[[1]], parameterised_wiring_terms_plot_list[[2]],
  pgs_distribution_plot, common.legend = FALSE, legend = "none", ncol = 3, align = "hv")
ggsave(filename = 'parameterised_gene_ontology_func_enrichment_plot.png', height = 5, 
       width = 12, dpi = 700)

### Part 6 - Chord diagram of seed network! ####
# Get the Schaefer 100-node regional labels
schaefer100_labels = unlist(schaefer100_coords[[1]])
rownames(seed) = schaefer100_labels
colnames(seed) = schaefer100_labels
# Format the data frame for plotting as a chord diagram
seed_df = data.frame(value = as.vector(seed), stringsAsFactors = FALSE,
  from = rep(rownames(seed), times = ncol(seed)),
  to = rep(colnames(seed), each = nrow(seed))) %>%
  # Remove any connections which are missing
  filter(value == 1)

chordDiagram(seed_df)

seed_df = data.frame(seed)


# This script visualizes the methodology schematic! Written by Alicja Monaghan
# in January 2025.

rm(list=ls())
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
library(cowplot)
library(grid)
library(rhdf5)
library(gridExtra)

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
eta_gamma_participant_level = read.csv('data/optimal_participant_level_gnm_parameters.txt')
# Load the nodal AHBA expression data
ahba_expression = read.csv('data/AHBA_expression_rnaseq_schaefer100_cleaned.csv')
# Load the functional enrichment results for parameterized nodal characteristics
nodal_costs_func_enrichment = read.csv(
  'data/nodal_wiring_costs_gost_formatted_results.csv')
# And the equivalent for parameterized nodal wiring value
nodal_value_func_enrichment = read.csv(
  'data/nodal_wiring_value_gost_formatted_results.csv'
)
# Load the group energy landscape
group_energy = data.frame(readMat('data/group_gnm_energy_99856.mat')[['original.simulation.energy']])
colnames(group_energy) = c("sptl", "neighbors", "matching", "clu-avg", "deg-avg")
# And get the parameters tested
group_parameters = read.csv('data/group_gnm_params_grid.csv')
# Load the densities of the thresholded target connectomes
schaefer100_thresholded_density = h5read(
  "data/ABCD_Individual_Target_Connectomes_Densities_schaefer100x17_Parcellation.mat",
  "ABCD_Thresholded_27_Streamlines_Density") %>% as.data.frame() %>%
  # Remove participant 1145 whose connectome reconstruction was incorrect (with
  # density of .2%).
  filter(V1 > 1)

### Part 1 - Visualizing seed and consensus network ####
png(filename = "seed_and_consensus_visualisation.png", width=5, height=5, units="in", res=700)
seed_graph = graph_from_adjacency_matrix(seed+consensus, 'undirected', weighted=TRUE)
E(seed_graph)$color[E(seed_graph)$weight > 1] = "darkgreen"
E(seed_graph)$color[E(seed_graph)$weight == 1] = "grey"
plot(seed_graph, vertex.size=0, vertex.shape='none', layout=data.matrix(coords),
     margin=-1, vertex.label = NA, axes = FALSE, label.cex=20)
dev.off()
# Find the density of each matrix!
print(paste("Density of seed matrix is", (sum(seed)/length(seed))*100, "%"))
print(paste("Density of consensus matrix is", (sum(consensus)/length(consensus))*100, "%"))
# Visualize the distribution of densities of thresholded connectomes
thresholded_connectome_density_boxplot = 
  ggplot(data=schaefer100_thresholded_density, mapping=aes(x = V1)) +
  geom_boxplot(outlier.shape=3, colour="#2F0F3E", width=.95) + labs(
    x=expression("Thresholded connectome density "*rho*" % (N = 2153)")) +
  theme(panel.background = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.x = element_line(color="black"),
        axis.text.x = element_text(size=45), axis.title.x = element_text(size=30),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1)) +
  scale_x_continuous(breaks = seq(2, 8, 2), limits = c(2, 8)) +
  scale_y_continuous(limits=c(-1, 1))
ggsave('thresholded_connectome_density_boxplot.png', plot = thresholded_connectome_density_boxplot,
       dpi = 700, height = 5, width = 10)

### Part 2 - Energy landscape for neighbors model ####
# Select the neighbors energy landscape, and merge with the parameters
neighbours_energy_df = data.frame(
  energy = group_energy$neighbors, eta = group_parameters$eta, 
  gamma = group_parameters$gamma)
# Specify the breaks for the x and y axes
axis_breaks = c(seq(from=-7, to=7, by=3.5))
neighbors_energy_landscape = 
  ggplot(neighbours_energy_df, aes(eta, gamma, fill = energy)) + geom_raster(interpolate=TRUE) +
  scale_fill_distiller(palette = "Spectral", name = "Energy") +
  labs(x = expression(eta), y = expression(gamma), title = "Homophily group-level energy",
       subtitle = "N = 99,856 runs") + 
  scale_x_continuous(breaks = axis_breaks) +
  # Flip the Y axis so that negative values are at the top
  scale_y_reverse(breaks=axis_breaks) +
  theme_classic() + theme(axis.line = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size=30), plot.title = element_text(size=27, hjust=0.5),
    axis.ticks = element_blank(), axis.title = element_text(size=35), text = element_text(size=35),
    legend.position = "left") +
  annotate("rect", xmin = -4.5, xmax = -1, ymin = 0, ymax = 0.5, alpha=.1, color="black") +
  coord_fixed() +
  guides(fill=guide_colorbar(ticks.colour=NA))
ggsave('neighbours_energy_landscape.png', dpi=700, height=6, width=8, plot=neighbors_energy_landscape)

### Part 3 - Distribution of optimal eta and gamma parameters (N = 2153) ####
eta_gamma_distribution_plot_main = ggplot(
  eta_gamma_participant_level, aes(x = optimal_eta, y =optimal_gamma)) +
  geom_pointdensity() + scale_color_paletteer_c("pals::ocean.matter") +
  labs(x = expression(eta), y = expression(gamma), title = "Optimal participant-level parameters",
       subtitle = "N = 74,529 runs") +
  theme(legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(color="black"), plot.title = element_text(size=25, hjust=0.5),
        text = element_text(size=35), axis.title = element_text(size=35),
        plot.subtitle = element_text(hjust=0.5, size=30))
# Add box plots showing the distribution of eta and gamma on the axes!
xdens = axis_canvas(eta_gamma_distribution_plot_main, axis = "x") +
  geom_boxplot(data=eta_gamma_participant_level, aes(x = optimal_eta),
               outlier.shape = 3, outlier.size = 2, color = "#2F0F3E")
ydens = axis_canvas(eta_gamma_distribution_plot_main, axis="y", coord_flip = TRUE) +
  geom_boxplot(data=eta_gamma_participant_level, aes(x = optimal_gamma),
               outlier.shape = 3, outlier.size = 2, color = "#2F0F3E") + coord_flip()
p1 = insert_xaxis_grob(eta_gamma_distribution_plot_main, xdens, grid::unit(.2, "null"), position = "top")
p2 = insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position="right")
# Combine with the neighbors energy landscape plot
ggsave('individual_parameter_estimates.png', dpi=700, width=8, height=6, plot=p2)

png(filename = 'group_energy_and_individual_parameter_estimates.png', width=10, height=5, units="in", res=700)
grid.arrange(grobs=list(neighbors_energy_landscape, p2), ncol=2)
dev.off()

### Part 4 - Polygenic scores for general cognitive ability ####
# Multiply the polygenic scores by 10000 for plotting the density
abcd_pgs_df$SCORE = abcd_pgs_df$SCORE*10000
pgs_distribution_plot = ggplot(
  data=abcd_pgs_df, mapping=aes(x = SCORE, y = PHENO)) + geom_pointdensity() + 
  labs(x = bquote("Polygenic scores for cognitive ability "(x10^-4)),
       y = "Cognitive ability") +
  theme(panel.background = element_blank(), legend.position = "none",
        axis.text = element_text(size=35), axis.title = element_text(size = 30),
        axis.line = element_line(color="black")) +
  scale_color_paletteer_c("pals::ocean.matter")
ggsave('abcd_intelligence_pgs_distribution.png', dpi=700, height=10, width=10)

### Part 5 - Visualize nodal distribution of AHBA regions previously linked to structural brain development ####
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
### Part 6 - Visualize gene enrichment of nodal wiring parameters ####
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
          axis.title = ggtext::element_markdown(size=25), axis.text = element_text(size=25),
          plot.title = element_text(size=25, hjust=0.5), legend.position = "none") +
    scale_colour_manual(values = c("GO:BP" = "#2F0F3E", "GO:CC" = "#CE4356", "GO:MF" = "#F29567")) +
    # Label the most important hits i.e. those with the smallest 
    geom_text_repel(data = subset(wiring_df, 0.25 < GeneRatio & GeneRatio < 1), aes(label=Description),
                    box.padding = 0.5, point.padding = 0.5, min.segment.length = 0, size = 8)
}
# Combine the two plots!
parameterised_gene_ontology_func_enrichment_plot = ggarrange(
  parameterised_wiring_terms_plot_list[[1]],
  parameterised_wiring_terms_plot_list[[2]], common.legend = FALSE, 
  legend = "none", ncol = 2, align = "hv")
ggsave(filename = 'parameterised_gene_ontology_func_enrichment_plot.png', height = 8, 
       width = 15, dpi = 700)

### Part 7 - Visualize wiring distance penalties by PGS group ####
# Identify the top and bottom 10% PGS groups
nsub = nrow(gnm_pgs_df)
bottom_pgs_boundary = sort(gnm_pgs_df$pgs)[round(nsub*.10)]
top_pgs_boundary = sort(gnm_pgs_df$pgs)[nsub-round(nsub*.10)]
gnm_pgs_df$pgs_level = 
  as.factor(ifelse(gnm_pgs_df$pgs<=bottom_pgs_boundary,"Lowest",
                   ifelse(gnm_pgs_df$pgs>=top_pgs_boundary, "Highest","NA")))
# Subset by those in either the top or bottom 10% of the PGS distribution
gnm_pgs_df_subsets = gnm_pgs_df %>% filter(str_detect(pgs_level, 'Lowest|Highest'))
# Plot the distribution of the wiring distance penalty eta by group
gnm_pgs_subset_eta_plot = ggplot(
  data=gnm_pgs_df_subsets, mapping=aes(x = eta, group=pgs_level, colour=pgs_level)) + 
  geom_boxplot() + coord_flip() + scale_colour_manual(values = c("#440154FF", "#1F968BFF")) +
  labs(x = expression(eta)) + scale_x_continuous(breaks = seq(-3.85, -2.25, by=.40)) +
  theme(legend.position = "none", panel.background = element_blank(),axis.ticks.x = element_blank(),
        axis.line.y = element_line(), axis.text.x = element_blank(), axis.text.y = element_text(size=25),
    axis.title.y = element_text(size=30))
ggsave('gnm_pgs_subset_eta_plot.png', dpi=700, width=8, height=8)

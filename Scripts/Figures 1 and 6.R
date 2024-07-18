# This script visualizes the relationship between polygenic scores for cognitive
# ability and generative network modelling parameters, alongside visualising 
# local and nodal properties of thresholded structural connectomes in ABCD.

### STEP 1 - Setting Up the Work Space ####
# Clear the workspace and set the working directory
rm(list = ls())
Working_Directory = "//cbsu/data/Imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_variation_structural_generative_mechanisms_open/"
setwd(Working_Directory)

# Now import the necessary packages
library(ggplot2)
library(reshape)
library(ggpointdensity)
library(ggpubr)
library(R.matlab)
library(tidyverse)
library(ggseg)
library(ggsegSchaefer)
library(matrixStats)
library(reshape2)
library(moments)
library(R.matlab)
library(rhdf5)
# Load the graph theory metrics of empirical connectomes
individual_gnm_targets_graph_theory_metrics = readMat("data/individual_gnm_targets_graph_theory_metrics.mat")
individual_gnm_targets_graph_theory_metrics = individual_gnm_targets_graph_theory_metrics[["abcd.individual.connectome.properties"]]
# And load the Schaefer 100-node parcellation labels.
schaefer100_labels = readMat("data/schaefer100x17_1mm_info.mat")
schaefer100_labels = unlist(schaefer100_labels[["schaefer100x17.1mm.info"]][[1]])

### PART 2 - Visualize distribution of local graph theory metrics ####
# Calculate mean nodal graph theory metrics across participants
mean_local_metric_tibble = tibble(
  degree = colMeans(individual_gnm_targets_graph_theory_metrics[,,1]),
  edge.length = colMeans(individual_gnm_targets_graph_theory_metrics[,,2]),
  betweenness.centrality = colMeans(individual_gnm_targets_graph_theory_metrics[,,3]),
  clustering = colMeans(individual_gnm_targets_graph_theory_metrics[,,4])
)
# Rescale all variables to have a range between 0 and 1
mean_local_metric_tibble = data.frame(
  lapply(mean_local_metric_tibble, function(x) 
    scale(x, center=FALSE, scale=max(x, na.rm = TRUE)/1)))
# Save the local graph theory metrics so that we can test their relationship 
# with the canonical sensorimotor-association axis developed by Sydnor and 
# colleagues (2021).
write.csv(x = mean_local_metric_tibble,file = 'mean_local_metric_tibble.csv')
# Reshape so that we have one column for the measure name and another for values
mean_local_metric_tibble = melt(mean_local_metric_tibble)
# Convert variable column to a factor
mean_local_metric_tibble$variable = factor(mean_local_metric_tibble$variable)
# Plot the distribution of metrics using box plots...
local_metrics_boxplot = 
  ggplot(mean_local_metric_tibble, aes(x = value, y = factor(variable), colour = factor(variable))) +
  geom_boxplot(outlier.colour = NULL, outlier.shape = 3, outlier.size = 2, linewidth = 1) +
  labs(x = 'Scaled Nodal Metric')+
  theme(legend.position = "none", panel.background = element_blank(),
        axis.line.x = element_line(colour="black"), axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20), axis.ticks.length.y = unit(0, "cm"),
        axis.title.x = element_text(size = 20), axis.ticks.length.x = unit(.15, "cm"),
        axis.title.y = element_blank()) +
  scale_color_manual(values = c('#440154FF', '#414487FF','#2A788EFF','#22A884FF')) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm"))
ggsave(filename = 'local.metrics.distribution.png', plot = local_metrics_boxplot,
       height = 5, width = 5, dpi = 700)
# Report the skewness for each metric.
local_metrics = unique(mean_local_metric_tibble$variable)
for (idx in 1:length(local_metrics)){
  metric_values = mean_local_metric_tibble %>% 
    filter(variable == local_metrics[idx]) %>% select('value')
  skewness_value = skewness(metric_values)
  print(sprintf('Skewness for %s is %.2f', local_metrics[idx], round(skewness_value[1], digits = 2)))
}

### PART 3 - Visualize distribution of global network properties ####
# Load the global graph theory metrics
global_metrics = read.csv('data/global_graph_theory_metrics.csv')
# Re-scale each column to have a range between 0 and 1
global_metrics = data.frame(
  lapply(global_metrics, function(x) 
    scale(x, center=FALSE, scale=max(x, na.rm = TRUE)/1)))
# Melt and convert variable to a factor!
global_metrics = melt(global_metrics) %>% mutate(variable = factor(variable))
# And visualize distribution using box plots
global_metrics_boxplot = 
  ggplot(data = global_metrics, mapping = aes(x = value, color = variable)) +
  geom_boxplot(outlier.colour = NULL, outlier.shape = 3, outlier.size = 2, linewidth = 1) +
  labs(x = 'Scaled Global Metric') +
  theme(legend.position = "none", panel.background = element_blank(),
        axis.line.x = element_line(colour="black"), axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20), axis.ticks.length.y = unit(0, "cm"),
        axis.title.x = element_text(size = 20)) +
  scale_color_manual(values = c('#440154FF', '#414487FF','#2A788EFF','#22A884FF')) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm"))
ggsave(filename = 'global.metrics.distribution.png', plot = global_metrics_boxplot,
       height = 5, width = 5, dpi = 700)
# Report the skewness for each metric.
global_metric_names = levels(global_metrics$variable)
for (idx in 1:length(global_metric_names)){
  metric_values = global_metrics %>% 
    filter(variable == global_metric_names[idx]) %>% pull('value')
  skewness_value = skewness(metric_values)
  print(sprintf('Skewness for %s is %.2f', global_metric_names[idx], round(skewness_value[1], digits = 2)))
}

### PART 4 - Visualize nodal coefficient of variation ####
cv_metric_array = array(NA, dim = c(100, length(local_metrics)))
for (local_metric_idx in 1:length(local_metrics)){
  # Calculate the metric's mean and standard deviation, across participants
  metric_mean = colMeans(individual_gnm_targets_graph_theory_metrics[,,local_metric_idx])
  metric_sd = colSds(individual_gnm_targets_graph_theory_metrics[,,local_metric_idx])
  metric_cv = metric_sd/metric_mean
  cv_metric_array[,local_metric_idx] = metric_cv
  # Create a data frame, and add the Schaefer 100-node 17-network labels
  variability_df = tibble(region = schaefer100_labels, cv = metric_cv)
  # Plot the distribution of variability across participants
  variability_brain_visualisation = 
    ggseg(.data = variability_df, atlas = schaefer17_100, 
        mapping = aes(fill = cv), hemisphere = "left") +
    theme_void() +
    labs(fill="") +
    paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1) +
    theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))
  ggsave(filename = sprintf('%s.variability.brain.plot.png', local_metrics[local_metric_idx]), 
         height = 2.5, width = 2.5, dpi=700, plot = variability_brain_visualisation)
}
cv_metric_df = data.frame(cv_metric_array)
colnames(cv_metric_df) = local_metrics
cv_metric_df$region = schaefer100_labels
### PART 5 - Visualize relationship between PGSs, GNM parameters, and Global Efficiency ####
# Visualise relationship between PGSs and eta...
eta_plot = ggplot(data = pgs_gnm, mapping=aes(x = pgs, y = eta)) +
  geom_pointdensity(size=2) +
  geom_smooth(method='lm', colour="black", se = FALSE) +
  labs(x = bquote("Polygenic Scores for Cognitive Ability "(x10^-4)),
       y = expression("Distance Penalty "~eta)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size=15),
        panel.background = element_blank(), axis.line = element_line(colour="black"),
        legend.position="bottom") +
  paletteer::scale_color_paletteer_c("pals::ocean.matter") +
  guides(color=guide_colorbar(ticks.colour = NA), title = "Density")
# Visualise relationship between PGSs and gamma...
gamma_plot = ggplot(data = pgs_gnm, mapping=aes(x = pgs, y = gamma)) +
  geom_pointdensity(size=2) +
  geom_smooth(method='lm', colour="black", se = FALSE) +
  labs(x = bquote("Polygenic Scores for Cognitive Ability "(x10^-4)),
       y = expression("Wiring Value "~gamma)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size=15),
        panel.background = element_blank(), axis.line = element_line(colour="black"),
        legend.position="bottom") +
  paletteer::scale_color_paletteer_c("pals::ocean.matter") +
  guides(color=guide_colorbar(ticks.colour = NA))

# Combine both plots and save
combined_plots = ggarrange(eta_plot, gamma_plot, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave(filename = 'combined_eta_gamma_relationship_with_pgs.png', dpi = 700, height = 5, width = 10)

# This script visualizes the relationship between polygenic scores for cognitive
# ability and generative network modelling parameters, alongside the
# relationship between generative parameters and network properties. We also
# visualize local and nodal properties of thresholded structural connectomes in
# ABCD, properties of simulated connectomes, and how generative network model 
# parameters vary across the population. Written by Alicja Monaghan in July 
# 2024, and updated in February 2025. 

### STEP 1 - Setting Up the Work Space ####
# Clear the work space and set the working directory
rm(list = ls())
setwd('/Users/alicjamonaghan/Desktop/abcd_genes_brain_cognition/')

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
library(rhdf5)
library(pals)
library(magrittr)

setwd('/Users/alicjamonaghan/Desktop/abcd_genes_brain_cognition/')
# Load the graph theory metrics of empirical connectomes
individual_gnm_targets_graph_theory_metrics = readMat("data/individual_gnm_targets_graph_theory_metrics.mat")
individual_gnm_targets_graph_theory_metrics = individual_gnm_targets_graph_theory_metrics[["abcd.individual.connectome.properties"]]
# And load the Schaefer 100-node parcellation labels.
schaefer100_labels = readMat("data/schaefer100x17_1mm_info.mat")
schaefer100_labels = unlist(schaefer100_labels[["schaefer100x17.1mm.info"]][[1]])
# Load the global graph theory metrics for the subset with PGSs
pgs_global_graph_theory_metrics = read.csv('data/global_graph_theory_metrics.csv')
# Load the unstandardized GNM parameters and co-variates for the 1399 
# participants, with eta estimates within 2 standard deviations of the mean. 
pgs_gnm = read.csv('data/filtered_pgs_gnm_unstandardized.csv')
# Load the simulated and observed connectome properties (lowest energy simulation)
across_participants_statistics = read.csv("data/across_participants_statistics.csv",sep = ",",header = F)
nodal_statistics_names = c(
  "Degree", "Clustering", "Betweenness centrality", "Edge length",
  "Local efficiency", "Eigenvector centrality", "Participation coefficient")
# Load the mean parameterised nodal wiring costs and values, respectively, 
# computed at an individual-level.
mean_parameterised_nodal_wiring_costs = c(
  h5read("data/Mean_Neighbour_Parameterised_Nodal_Cost_Across_Participants.mat",
        "/Mean_Neighbour_Parameterised_Nodal_Cost"))
mean_parameterised_nodal_wiring_value = c(
  h5read("data/Mean_Neighbour_Parameterised_Nodal_Value_Across_Participants.mat",
         "/Mean_Neighbour_Parameterised_Nodal_Value"))
# Load the energy values for each participant's best simulation
individual_energy = c(h5read("data/ABCD_Individual_Lowest_Energy.mat", "/ABCD_Individual_Lowest_Energy"))
# And the optimal individual-level eta and gamma estimates
individual_gnm_parameters = read.table("data/optimal_participant_level_gnm_parameters.txt", header = TRUE, sep = ",")
# Load the global efficiency for the attenuated network parameters
efficiency_attenuation = readMat(
  'data/synthetic_connectome_decile_global_efficiency.mat')[["synthetic.connectome.global.efficiency"]] %>%
   t() %>% as.data.frame() %>% `colnames<-` (
     c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90")) %>%
  melt()
# Load the topological dissimilarity for simulated networks versus randomly-
# rewired networks across attenuated network parameters
topological_dissimilarity_attenuation = readMat(
  'data/topological_dissimilarity_random_rewiring.mat')[["rand.dissim.pct"]] %>%
  t() %>% as.data.frame() %>% `colnames<-` (
    c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90")) %>% melt()
# Load the data we have pulled for Duncan's exploratory analysis
data_for_duncan = read.csv('data/data_for_duncan.csv')

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
# Remove the participant IDs when melting the global graph theory metrics data frame
global_metrics = subset(pgs_global_graph_theory_metrics, select=-c(id))
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
### PART 5 - Visualize relationship between PGS and GNM parameters ####
# Multiply PGS scores by 10000 for visualisation purposes
pgs_gnm$pgs = pgs_gnm$pgs * 10000
# Visualize relationship between PGSs and eta...
eta_plot = ggplot(data = pgs_gnm, mapping=aes(x = pgs, y = eta)) +
  geom_pointdensity(size=2) +
  geom_smooth(method='lm', colour="darkgrey", se = FALSE) +
  labs(x = bquote("Polygenic scores for cognitive ability "(x10^-4)),
       y = expression("Distance penalty "~eta)) +
       #title = "Softer wiring distance penalties are associated \nwith stronger polygenic scores for cognitive ability")+
   theme(axis.text = element_text(size = 15), axis.title = element_text(size=15),
        panel.background = element_blank(), axis.line = element_line(colour="black"),
        legend.position="none", plot.title = element_text(size=15, hjust=0.5)) +
  paletteer::scale_color_paletteer_c("pals::ocean.matter") 
# Visualize relationship between PGSs and gamma...
gamma_plot = ggplot(data = pgs_gnm, mapping=aes(x = pgs, y = gamma)) +
  geom_pointdensity(size=2) +
  geom_smooth(method='lm', colour="darkgrey", se = FALSE) +
  labs(x = bquote("Polygenic scores for cognitive ability "(x10^-4)),
       y = expression("Wiring value "~gamma))+
       #title = "Wiring value is not associated with \npolygenic scores for cognitive ability") +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size=15),
        panel.background = element_blank(), axis.line = element_line(colour="black"),
        legend.position="none", plot.title = element_text(size=15, hjust=0.5)) +
  paletteer::scale_color_paletteer_c("pals::ocean.matter") 

### PART 6 - Correlations between observed and simulated nodal properties ####
# Create an array to hold the Pearson correlation coefficients, alongside the 
# uncorrected p-values and Bonferroni-corrected p-values
correlation_results_array = array(NA, dim=c(length(nodal_statistics_names), 3))
# Correlate the observed and simulated connectome properties
for (i in 1:length(nodal_statistics_names)){
  correlation_result = cor.test(
    across_participants_statistics[, i], across_participants_statistics[, i+8],
    method="pearson")
  # Assign to the output array
  correlation_results_array[i, 1] = correlation_result$estimate
  correlation_results_array[i, 2] = correlation_result$p.value
}
# Correct for multiple comparisons using the Bonferroni method
correlation_results_array[, 3] = p.adjust(correlation_results_array[, 2], method="bonferroni")
# Convert the results array to a data frame and categorize the adjusted p-values
# according to whether they are below .05
correlation_results_df = correlation_results_array %>% as.data.frame() %>% 
  setNames(c("coefficient", "unadjusted_p", "adjusted_p")) %>%
  mutate(measure = nodal_statistics_names) %>%
  # Take the absolute correlation coefficient - the only slight negative 
  # coefficient is betweenness-centrality: we'll note this in the legend.
  mutate(coefficient = abs(coefficient)) %>%
  mutate(significance_category = factor(
    ifelse(adjusted_p <= .05, "Significant", "Non-Significant")))
# Plot the correlation coefficients as circles: these are filled if the adjusted
# p-value is below .05. Order by decreasing correlation coefficient.
correlation_scatter_plot = ggplot(
  data=correlation_results_df, mapping=aes(x=reorder(measure, coefficient), y=coefficient, color=significance_category, fill=significance_category)) +
  geom_bar(stat="identity", alpha=.6, width=.5) + 
  labs(x = "Nodal measure", y = "Pearson correlation coefficient between simulated\n and observed connectomes",
       title = "Participant-level topological features captured by\n lowest energy simulations",
       subtitle = "N = 2153") +
  coord_flip() + theme(
    panel.background = element_blank(),axis.line = element_line(), axis.ticks.y = element_blank(),
    axis.text = element_text(size=15), axis.title = element_text(size=15), legend.position = "none",
    plot.title = element_text(size=15, hjust=0.5), plot.subtitle = element_text(size=15, hjust=0.5),
    axis.ticks.length.x = unit(1, "mm"), plot.margin = margin(.5,.5,.5,.5, "cm")) +
  scale_y_continuous(expand=c(0,0), limits = c(0, .8), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8)) +
  scale_color_manual(values = c("Non-Significant" = "lightgrey", "Significant" = "#9D2461")) +
  scale_fill_manual(values = c("Non-Significant" = "lightgrey", "Significant" = "#9D2461"))
ggsave(correlation_scatter_plot, filename = 'correlation_scatter_plot.png', dpi=700, height = 5, width = 8)
### PART 7 - Plot parameterized wiring costs and values averaged across participants ####
# Create an array with the parameterised nodal wiring costs, values, and labels
parameterised_df = data.frame(
  costs = mean_parameterised_nodal_wiring_costs,
  values = mean_parameterised_nodal_wiring_value,
  region = schaefer100_labels)
parameterised_names_to_plot = c("costs", "values")
# Plot the parameterized nodal wiring costs and values across the left hemisphere.
for (i in 1:length(parameterised_names_to_plot)){
  parameterised_brain_visualisation = 
    ggseg(.data = parameterised_df, atlas = schaefer17_100, hemisphere="left",
        mapping=aes(fill = .data[[parameterised_names_to_plot[i]]])) +
    theme_void() + theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1)
  ggsave(filename = sprintf('parameterised_nodal_wiring_%s_brain_plot.png', parameterised_names_to_plot[i]),
         height = 2.5, width = 2.5, dpi=700, plot = parameterised_brain_visualisation)
}
### PART 8 - Plot distribution of global efficiency across attenuations ####
efficiency_attenuation_plot = ggplot(
  data=efficiency_attenuation, mapping=aes(x=variable, y=value, color=variable)) +
  geom_boxplot(outlier.colour = NULL, outlier.shape = 3, outlier.size = 2) +
  labs(x = "Percentage of optimal parameters (%)", y = "Global efficiency") +
       #title = "Relationship between stochasticity and global efficiency") +
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=15), axis.title = element_text(size=15),
        legend.position = "none", plot.title = element_text(size=15, hjust=0.5)) +
  # Reverse order of X axis i.e. start at 90% and end at 0%
  scale_x_discrete(limits=rev(levels(efficiency_attenuation$variable))) +
  scale_y_continuous(limits=c(0.25, .45), n.breaks = 5) +
  # Use the same color scheme as the other plots
  scale_colour_manual(values=rev(unname(ocean.matter(n=10))))

### PART 9 - Plot distribution of topological dissimilarity to random networks ####
TD_attenuation_plot = ggplot(
  data=topological_dissimilarity_attenuation, mapping=aes(x=variable, y=value, color=variable)) +
  geom_boxplot(outlier.colour=NULL, outlier.shape=3, outlier.size=2) +
  labs(x = "Percentage of optimal parameters (%)", y = expression(Delta~"Topological fingerprint"))+
      #title = "Dissimilarity to random networks") +
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=15), axis.title = element_text(size=15),
        legend.position = "none", plot.title = element_text(size=15, hjust=0.5)) +
  # Reverse order of X axis i.e. start at 90% and end at 0%
  scale_x_discrete(limits=rev(levels(topological_dissimilarity_attenuation$variable))) +
  scale_colour_manual(values=rev(unname(ocean.matter(n=10))))

### PART 10 - Compile plots ####
compiled_plot = ggarrange(
  plotlist = list(eta_plot, gamma_plot, efficiency_attenuation_plot,
                  TD_attenuation_plot), ncol = 4, align="h")
ggsave(filename='eta_gamma_simulations_global_efficiency.png', compiled_plot,
       height = 6, width=25, dpi=700)


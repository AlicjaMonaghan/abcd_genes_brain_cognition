# This code creates key figures displaying the relationship between generative
# network modelling (GNM) parameters, cognitive ability polygenic scores (PGSs),
# region-wise gene expression values from the Allen Human Brain Atlas (AHBA),
# and a general g factor of cognitive ability, from 2153 children in the 
# baseline time point of the fourth release of the Adolescent Brain Cognitive 
# Development (ABCD) Study. Correspondence to Alicja Monaghan, 
# alicja.monaghan@mrc-cbu.cam.ac.uk.

# STEPS:
# 1. Setting up the work space.
# 2. Loading the data for the primary manuscript visualisation.
# 3. Plotting nodal graph theory properties (mean degree, total edge length, 
# betweenness-centrality, and clustering) across participants for GNM targets 
# (FIGURES 1D-G). These targets were individual structural connectomes 
# thresholded at 27 streamlines.
# 4. Plotting energy distributions for the group-level GNMs in the Schaefer 100-
# node (17-network) parcellation (FIGURE 2D).
# 5. Plotting distribution of nodal degree for simulated and observed 
# connectomes (FIGURE 2F). Note that we visualise simulated degree 
# representative of the mean correlation between simulated and observed degree 
# across 1000 simulations of the lowest-energy eta-gamma combination for each 
# rule.
# 6. Plotting mean parameterised nodal wiring costs and values from individual
# GNMs for 2153 children (FIGURES 3B-C).
# 7A. Visualising relationship between simulated and observed local statistics
# from individual-level GNMs (FIGURES 3D-K) e.g. participation coefficient, 
# and optimal community structure given modularity.
# 7B. Plotting regional largest and smallest errors in simulated local 
# statistics, showing over/under-stepping. (FIGURES 3D-K).
# 8. Chord diagram showing bootstrapped X loadings of a PLS with eta, gamma,
# and PGSs as predictors of cognitive ability across 1461 children (FIGURE 4A).
# 9A. Loading data of the mean number of samples provided by each AHBA donor for
# each region across 3 parcellations.
# 9B. Plotting these mean distributions - SUPPLEMENTARY FIGURE 1.

### PART 1 - Set Up the Work Space ####
rm(list = ls())
library(R.matlab)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggseg)
library(ggsegSchaefer)
library(raveio)
library(viridis)
library(circlize)
library(ggseg3d)
library(ggsegBrainnetome)
# SET WORKING DIRECORY TO WHERE YOU SAVED THIS DIRECTORY!
setwd('abcd_genomic_variation_structural_generative_mechanisms_open/')

### PART 2 - Load Data for Main Visualisation ####
# Load the group-level GNM energy for 5 wiring rules and 99,856 eta-gamma combinations.
group_gnm_energy = readMat("data/group_gnm_energy_99856.mat")
# Unlist. 
group_gnm_energy = data.frame(group_gnm_energy[["original.simulation.energy"]])
# Set the order of models in the array. 
modeltypes = c("Spatial", "Neighbours", "Matching", "Clu-Avg", "Deg-Avg")
# And format column names.
colnames(group_gnm_energy) = modeltypes
# Set the number of group simulations we used.
nruns = nrow(group_gnm_energy)
runs_to_visualise = 1000
# Load the nodal degree distributions for observed and simulated networks. Note,
# the simulated degrees are averaged across 1000 runs of the lowest-energy 
# parameter combination for each model.
group_gnm_degree_distributions = readMat("data/degree_simulated_group_gnms.mat")
group_gnm_degree_distributions = data.frame(t(group_gnm_degree_distributions[["degree.to.plot"]]))
# Load nodal degree, edge length, betweenness-centrality, and clustering 
# coefficient values for individual GNM targets.
individual_gnm_targets_graph_theory_metrics = readMat("data/individual_gnm_targets_graph_theory_metrics.mat")
individual_gnm_targets_graph_theory_metrics = individual_gnm_targets_graph_theory_metrics[["abcd.individual.connectome.properties"]]
# And load the Schaefer 100-node parcellation labels.
schaefer100_labels = readMat("data/schaefer100x17_1mm_info.mat")
schaefer100_labels = unlist(schaefer100_labels[["schaefer100x17.1mm.info"]][[1]])
# Load the mean parameterised nodal wiring costs and values, respectively, 
# computed at an individual-level.
mean_parameterised_nodal_wiring_costs = unlist(read_mat("data/Mean_Neighbour_Parameterised_Nodal_Value_Across_Participants.mat"))
mean_parameterised_nodal_wiring_value = unlist(read_mat("data/Mean_Neighbour_Parameterised_Nodal_Cost_Across_Participants.mat"))
# Load the simulated and observed statistics averaged across participants for 
# the individual-level GNMs. The first 8 columns represents the observed 
# statistics, whilst the final 8 columns represents corresponding simulated 
# statistics.
across_participants_statistics = read.csv("data/across_participants_statistics.csv",sep = ",",header = F)
# Load the summary table of the PLS with eta, gamma, and PGSs for cognitive 
# ability predicting cognitive ability. These are the bootstrapped X loadings 
# corrected for sex, in-scanner motion, site, and age as co-variates.
cognitive_ability_onto_eta_gamma_pgs_plus_covariates_bootstrapped_xloadings = 
  unlist(read_mat("data/cognitive_ability_onto_eta_gamma_pgs_plus_covariates_bootstrapped_xloadings.mat"))
### PART 3 - Plotting Nodal Graph Theory Properties Across Participants for GNM Targets - FIGURE 1D-G ####
# Plot the nodal distributions of 4 graph theory metrics for individual GNM
# targets i.e. connectomes thresholded at a 60% consensus, and then at 27 
# streamlines. Specify these metrics.
local_metrics_to_plot = 
  c("Mean Degree","Mean Edge Length","Mean Betweenness-Centrality","Mean Clustering coefficient")
# Specify the colour palette which we shall use. 
colour_palette = brewer.pal(n=9,"Reds")
# Create a list to hold the plots.
local_metric_plots_list = vector("list", length(local_metrics_to_plot))
# For each metric, average across participants to generate one value per node. 
# Create a tibble with the mean value and the labels for the Schaefer 100-node
# 17-network parcellation. Assign the plot to the output list.
for (local_metric in 1:length(local_metrics_to_plot)){
  # Create a tibble with the region labels and the nodal means. 
  mean_local_metric_tibble = tibble(region = schaefer100_labels,
                                    Mean = colMeans(individual_gnm_targets_graph_theory_metrics[,,local_metric]))
  # Plot!
  local_metric_plots_list[[local_metric]] = mean_local_metric_tibble %>%
    ggplot() +
    geom_brain(atlas = schaefer17_100, position = position_brain(hemi~side),
               aes(fill = Mean), colour = "darkgrey") +
    scale_fill_gradientn(colours = colour_palette) +
    ggtitle(local_metrics_to_plot[local_metric]) +
    theme_void() +
    theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
          axis.text = element_blank(), axis.ticks = element_blank(),
          legend.title = element_blank(), legend.key.width = unit(2.5,'cm'),
          legend.key.height = unit(1,'cm'),
          legend.text = element_text(size = 25),
          legend.position = "bottom", legend.direction = "horizontal") +
    scale_x_continuous(breaks = 6) +
    guides(fill = guide_colorbar(ticks.colour = NA))
  
}

### PART 4 - Plotting Energy Distribution - FIGURE 2D ####
# Extract lowest-energy 10% simulations from group_gnm_outputs and format.
group_gnm_energy_sorted = 
  tibble("Spatial" = sort(group_gnm_energy$Spatial)[1:runs_to_visualise],
         "Neighbours" = sort(group_gnm_energy$Neighbours)[1:runs_to_visualise],
         "Matching" = sort(group_gnm_energy$Matching)[1:runs_to_visualise],
          "Clu-Avg" = sort(group_gnm_energy$`Clu-Avg`)[1:runs_to_visualise],
          "Deg-Avg" = sort(group_gnm_energy$`Deg-Avg`)[1:runs_to_visualise])

# Melt the above data frame as we need a long-format for plotting.
group_gnm_energy_sorted = reshape2::melt(group_gnm_energy_sorted)
# Create a model class factor e.g. spatial, homophily, clustering, degree.
modelmembership = c("Spatial","Homophily","Clustering","Degree")
group_gnm_energy_sorted$'Model Membership' = 
  factor(rep(modelmembership, c(runs_to_visualise, runs_to_visualise*2, runs_to_visualise, runs_to_visualise)))
# Create a factor for each wiring rule and format column names.
colnames(group_gnm_energy_sorted)[1:2] = c("Model Type", "Energy")
group_gnm_energy_sorted$`Model Type` = as.factor(group_gnm_energy_sorted$`Model Type`)
# Now plot the energy distributions!
group_energy_distribution_plot = group_gnm_energy_sorted %>%
  ggplot(aes(x = `Model Type`,y = Energy, fill = `Model Membership`)) + 
  stat_boxplot(geom ='errorbar', width=.25) +
  geom_boxplot(width = .25, outlier.shape = NA) +
  labs(x = "Model Type", y = "Energy", 
       title = "Group Generative Modelling Energy Distributions",
       subtitle = paste0(runs_to_visualise," Lowest-Energy Simulations from N = ",nruns, " Simulations")) +
  theme(plot.title = element_text(size=30,hjust=0.5,face="bold"),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(face="bold", size = 25),
        legend.position = "top", 
        axis.text = element_text(size = 25), legend.text = element_text(size = 25),
        legend.key.size = unit(1,'cm'), legend.background = element_rect(fill='white', color = 'black'),
        legend.key = element_blank(),
        axis.ticks = element_blank()) +
  ylim(0,max(group_gnm_energy_sorted$Energy))
print(group_energy_distribution_plot)

### PART 5 - Plotting Nodal Degree Distributions of Simulations and Observed Connectomes - FIGURE 2F ####
# Whilst energy is one evaluative metric to choose the best GNM, where the lower
# the energy, regardless of distribution, the better the model. However, we also
# want to ensure that the GNM can capture local statistical properties, such as 
# the distribution of degree (number of connections each node has). Therefore, 
# we shall correlate nodal degree in the observed connectomes with each of the 
# 5 wiring rules. For each rule, we conducted 1000 simulations of the lowest-
# energy parametrr combination, and calculated the mean correlation between
# simulated and observed degree. We then selected the degree distribution 
# representative of this mean to plot here. 
# First, extract the observed and simulated nodal degree and format.
colnames(group_gnm_degree_distributions)[1:5] = modeltypes
colnames(group_gnm_degree_distributions)[6] = "Observed"

# For each of the simulated models, plot the Pearson correlation between their 
# nodal degree distribution with the observed connectome. 

# SPATIAL MODEL # 
cor.test(group_gnm_degree_distributions$Observed, group_gnm_degree_distributions$Spatial)
ggplot(group_gnm_degree_distributions, mapping = aes(x = Observed, y = Spatial)) +
  geom_point(size = 4, colour = "grey") +
  geom_smooth(method="lm", se = FALSE, color = "red", size = 2) +
  labs(x = "Observed Nodal Degree", y = "Simulated Nodal Degree", title = "Spatial",
       subtitle = expression(bold(bolditalic("r")~"(98) = .131, "~bolditalic("p")~" = .193"))) +
  theme(plot.title = element_text(size = 25, face="bold", hjust = .5),
        plot.subtitle = element_text(size = 25, face="bold", hjust = .5, colour = "red"),
        axis.title = element_text(size = 25, face = "bold", hjust = .5),
        axis.text = element_text(size = 30),
        panel.background = element_blank(), axis.line = element_line(colour="black"))
ggsave("Nodal_Degree_Observed_Spatial_Schaefer100.jpg", dpi = 700, height = 6, width = 6)

# HOMOPHILY-NEIGHBOURS # 
cor.test(group_gnm_degree_distributions$Observed, group_gnm_degree_distributions$Neighbours)
ggplot(group_gnm_degree_distributions, mapping = aes(x = Observed, y = Neighbours)) +
  geom_point(size = 4, colour = "grey") +
  geom_smooth(method="lm", se = FALSE, color = "red", size = 2) +
  labs(x = "Observed Nodal Degree", y = "Simulated Nodal Degree", title = "Neighbours",
       subtitle = expression(bold(bolditalic("r")~"(98) = .406, "~bolditalic("p")~" < .001"))) +
  theme(plot.title = element_text(size = 25, face="bold", hjust = .5),
        plot.subtitle = element_text(size = 25, face="bold", hjust = .5, colour = "red"),
        axis.title = element_text(size = 25, face = "bold", hjust = .5),
        axis.text = element_text(size = 30),
        panel.background = element_blank(), axis.line = element_line(colour="black"))
ggsave("Nodal_Degree_Observed_Neighbours_Schaefer100.jpg", dpi = 700, height = 6, width = 6)

# HOMOPHILY-MATCHING # 
cor.test(group_gnm_degree_distributions$Observed, group_gnm_degree_distributions$Matching)
ggplot(group_gnm_degree_distributions, mapping = aes(x = Observed, y = Matching)) +
  geom_point(size = 4, colour = "grey") +
  geom_smooth(method="lm", se = FALSE, color = "red", size = 2) +
  labs(x = "Observed Nodal Degree", y = "Simulated Nodal Degree", title = "Matching",
       subtitle = expression(bold(bolditalic("r")~"(98) = .393, "~bolditalic("p")~" < .001"))) +
  theme(plot.title = element_text(size = 25, face="bold", hjust = .5),
        plot.subtitle = element_text(size = 25, face="bold", hjust = .5, colour = "red"),
        axis.title = element_text(size = 25, face = "bold", hjust = .5),
        axis.text = element_text(size = 30),
        panel.background = element_blank(), axis.line = element_line(colour="black"))
ggsave("Nodal_Degree_Observed_Matching_Schaefer100.jpg", dpi = 700, height = 6, width = 6)

# AVERAGE-CLUSTERING # 
cor.test(group_gnm_degree_distributions$Observed, group_gnm_degree_distributions$`Clu-Avg`)
ggplot(group_gnm_degree_distributions, mapping = aes(x = Observed, y = `Clu-Avg`)) +
  geom_point(size = 4, colour = "grey") +
  geom_smooth(method="lm", se = FALSE, color = "red", size = 2) +
  labs(x = "Observed Nodal Degree", y = "Simulated Nodal Degree", title = "Clu-Avg",
       subtitle = expression(bold(bolditalic("r")~"(98) = .166, "~bolditalic("p")~" = .099"))) +
  theme(plot.title = element_text(size = 25, face="bold", hjust = .5),
        plot.subtitle = element_text(size = 25, face="bold", hjust = .5, colour = "red"),
        axis.title = element_text(size = 25, face = "bold", hjust = .5),
        axis.text = element_text(size = 30),
        panel.background = element_blank(), axis.line = element_line(colour="black"))
ggsave("Nodal_Degree_Observed_Clu_Avg_Schaefer100.jpg", dpi = 700, height = 6, width = 6)

# AVERAGE-DEGREE # 
cor.test(group_gnm_degree_distributions$Observed, group_gnm_degree_distributions$`Deg-Avg`)
ggplot(group_gnm_degree_distributions, mapping = aes(x = Observed, y = `Deg-Avg`)) +
  geom_point(size = 4, colour = "grey") +
  geom_smooth(method="lm", se = FALSE, color = "red", size = 2) +
  labs(x = "Observed Nodal Degree", y = "Simulated Nodal Degree", title = "Deg-Avg",
       subtitle = expression(bold(bolditalic("r")~"(98) = .399, "~bolditalic("p")~" < .001"))) +
  theme(plot.title = element_text(size = 25, face="bold", hjust = .5),
        plot.subtitle = element_text(size = 25, face="bold", hjust = .5, colour = "red"),
        axis.title = element_text(size = 25, face = "bold", hjust = .5),
        axis.text = element_text(size = 30),
        panel.background = element_blank(), axis.line = element_line(colour="black"))
ggsave("Nodal_Degree_Observed_Deg_Avg_Schaefer100.jpg", dpi = 700, height = 6, width = 6)

### PART 6 - Plot Mean Parameterized Nodal wiring Costs and Values - FIGURE 3B-C ####
# We shall plot the parameterised nodal wiring costs and values, averaged across
# participants to produce nodal data, on cortical surface using the Schaefer 
# 100-node parcellation. First, create a data frame compatible with ggseg.
mean_parameterised = tibble(region = schaefer100_labels,
                            nodal_costs = mean_parameterised_nodal_wiring_costs,
                            nodal_value = mean_parameterised_nodal_wiring_value)
# Now visualise the nodal costs
ggseg(mean_parameterised, atlas = schaefer17_100, colour = "black", size = 1, 
      mapping = aes(fill = nodal_costs), position = "stacked") +
  scale_fill_viridis(label = function(x) sprintf("%.3f", x)) +
  theme_void() +
  labs(title = "Parameterised Nodal Wiring Costs",
       fill = expression(sum(D["i,:"])^eta)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 20), 
        legend.text = element_text(size = 15)) +
  guides(fill=guide_colorbar(ticks.colour = NA))

# And visualise nodal value
ggseg(mean_parameterised, atlas = schaefer17_100, colour = "black", size = 1, 
      mapping = aes(fill = nodal_value), position = "stacked") +
  scale_fill_viridis(label = function(x) sprintf("%.3f", x)) +
  theme_void() +
  labs(title = "Parameterised Nodal Wiring Value",
       fill = expression(sum(K["i,:"])^gamma)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 20), 
        legend.text = element_text(size = 15)) +
  guides(fill=guide_colorbar(ticks.colour = NA))

### PART 7 - Visualize Relationship Between Simulated and Observed Local Statistics - FIGURE 3D-K ####
# This is made up of two sections. The first will plot the correlations between 
# simulated and observed local statistics. The second will plot the the largest
# and smallest regional errors in the simulations for these local statistics on
# a 3D rendering of the brain. 
### PART 7A - Plot Correlations Between Simulated and Observed Statistics ####
# First, create a vector of the nodal statistics we will plot.
nodal_statistics_names = 
  c("Degree (k)", "Clustering (c)", "Betweenness (b)", "Edge Length (d)", 
    "Local Efficiency (e)", "Eiegnvector Centrality (EC)", 
    "Participation coefficient (PC)", "Modularity (Q)")
# And the associated symbols alone
nodal_statistics_symbols = c("k", "c", "b", "d", "e", "EC", "PC", "Q")
# Initialize a list which will hold the observed vs simulated data frames for 
# each nodal statistic, as well as the error levels!
nodal_observed_simulated_dataframe_list = vector("list", length(nodal_statistics_names))
# Create a loop which will plot the simulated and observed nodal statistics. 
# Note that this loop will create data frames showing the simulated and observed 
# values, and the Schaefer 100-node labels, for each local statistic. 
for (variable_index in 1:length(nodal_statistics_names)){
  # Create a data frame showing the nodal simulated and observed statistics.
  nodal_observed_simulated_dataframe = 
    tibble(Observed = across_participants_statistics[,variable_index],
           Simulated = across_participants_statistics[,8+variable_index],
           region = schaefer100_labels,
           error = Observed-Simulated)
  # Calculate the correlation between the nodal observed and simulated statistics
  correlation = cor.test(nodal_observed_simulated_dataframe$Observed, nodal_observed_simulated_dataframe$Simulated, method = "pearson")
  # Format the correlation coefficient and p-value for plotting
  correlation_coefficient = round(correlation$estimate, 3)
  if (correlation$p.value < .001){
    correlation_P_Value = "< .001"
  } else {
    correlation_P_Value = paste("=", round(correlation$p.value, 3))
  }
  # Create a column which will indicate whether the regional error is the highest
  # (top 10%), lowest (bottom 10%), or in-between (NA).
  highest_error_boundary = sort(abs(nodal_observed_simulated_dataframe$error))[91]
  lowest_error_boundary = sort(abs(nodal_observed_simulated_dataframe$error))[10]
  nodal_observed_simulated_dataframe$'Error_Level' = 
    ifelse(abs(nodal_observed_simulated_dataframe$error) >= highest_error_boundary, "Largest",
           ifelse(abs(nodal_observed_simulated_dataframe$error) <= lowest_error_boundary, "Smallest", NA))
  # Assign the different error factor levels to different colours. These will be
  # used when plotting as a 3D mesh!
  nodal_observed_simulated_dataframe$'Error_Colour' = 
    ifelse(nodal_observed_simulated_dataframe$Error_Level=="Largest", "#C43C4EFF", 
           ifelse(nodal_observed_simulated_dataframe$Error_Level=="Smallest", "#56C667FF", NA))
  # Assign the fully formatted and complete data frame to the list.
  nodal_observed_simulated_dataframe_list[[variable_index]] = nodal_observed_simulated_dataframe
  # Plot the correlation between this simulated and observed nodal statistic. We
  # will colour the individual points according to their error level factor, so
  # that we can visualise to what extent the model is over- or under-estimating
  # each statistic.
  ggplot(nodal_observed_simulated_dataframe,
         mapping = aes(x = Observed, y = Simulated, colour = as.factor(Error_Level))) +
    geom_point(size = 4, show.legend = F) +
    geom_smooth(method="lm", se = FALSE, colour = "#453781FF", lwd = 2) + 
    labs(title = nodal_statistics_names[variable_index],
         subtitle = bquote(bolditalic("r")~bold("(98) = "~ .(as.character(correlation_coefficient)) ~ "," 
                                                ~ bolditalic("p") ~ .(as.character(correlation_P_Value)))),
         x = paste("Observed", nodal_statistics_symbols[variable_index]), 
         y = paste("Simulated", nodal_statistics_symbols[variable_index])) +
    theme(axis.title = element_text(size = 30, face = "bold", hjust = .5),
          axis.text = element_text(size = 30), axis.ticks = element_blank(),
          axis.line = element_line(colour = "black"), panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 30, face = "bold", hjust = .5),
          plot.subtitle = element_text(size = 30, hjust = .5, colour = "red")) +
    scale_color_manual(values = c("Largest" = "#C43C4EFF", "Smallest" = "#56C667FF", "NA" = "grey"))
  filename = paste0(nodal_statistics_names[variable_index],'.png')
  ggsave(filename, units="in", width=9, height=6, dpi = 700)
  # Remove unnecessary variables
  rm(nodal_observed_simulated_dataframe, correlation, correlation_coefficient, correlation_P_Value)
}
### PART 7B - Visualize Regions with Smallest and Largest Errors, Respectively ####
# Use the ggseg3d function to plot the lowest and largest regional errors for
# the left and right hemispheres, respectively. Note that some errors in the 
# ggseg3d function are raised, but this is mainly due to deprecated functions in
# tidyr 1.2.0.

### 3D DEGREE ###
ggseg3d(nodal_observed_simulated_dataframe_list[[1]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D CLUSTERING ###
ggseg3d(nodal_observed_simulated_dataframe_list[[2]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D BETWEENNESS CENTRALITY ###
ggseg3d(nodal_observed_simulated_dataframe_list[[3]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D EDGE LENGTH ###
ggseg3d(nodal_observed_simulated_dataframe_list[[4]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D LOCAL EFFICIENCY ###
ggseg3d(nodal_observed_simulated_dataframe_list[[5]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D EIGENVECTOR CENTRALITY ###
ggseg3d(nodal_observed_simulated_dataframe_list[[6]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D PARTICIPATION COEFFICIENT ###
ggseg3d(nodal_observed_simulated_dataframe_list[[7]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### 3D MODULARITY ###
ggseg3d(nodal_observed_simulated_dataframe_list[[8]], atlas = schaefer17_100_3d,
        surface = "inflated", label = "region", hemisphere = "left", 
        colour = "Error_Colour", show.legend = T, na.alpha = .5) %>% 
  remove_axes() %>% 
  pan_camera("left medial")

### PART 8 - Plot PLS with Eta, Gamma, and PGS as Predictors of Cognitive Ability - FIGURE 4A ####
# To explore to what extent linear combinations of structural brain development
# parameters (eta and gamma) and polygenic scores (PGS) for general cognitive 
# ability predicted linear variation in cognitive ability, we conducted a 
# partial least squares (PLS) regression with eta, gamma, and PGSs as predictors
# of cognitive ability. Note that we included sex, age, scanner site, and mean
# frame wise displacement as co-variates, and corrected for these through partial 
# correlation of X and Y scores. First, convert the bootstrapped X loadings into
# a matrix.
bootstrapped_xloadings = matrix(cognitive_ability_onto_eta_gamma_pgs_plus_covariates_bootstrapped_xloadings,3,3)
# Use unicode names for eta and gamma special characters.
rownames(bootstrapped_xloadings) = c("\U03B7","\U03B3", "PGS")
rownames(bootstrapped_xloadings) = c("\U1D6C8","\U1D6C4", "PGS")
colnames(bootstrapped_xloadings) = paste0("LV ",1:3)
# Visualise negative loadings in blue and positive in red, with white as the midpoint. 
col_fun = colorRamp2(c(min(bootstrapped_xloadings), 0, max(bootstrapped_xloadings)),
                     c("#377EB8","white","#E41A1C"), transparency = .5)
# Change the font size of the labels
par(cex = 2, mar = c(0, 0, 0, 0))
# Add gaps between the predictors and components, and rotate the plot so that 
# the predictors are on the left hand side. 
circos.par(gap.after = c(rep(5, nrow(bootstrapped_xloadings)-1), 15, 
                         rep(5, ncol(bootstrapped_xloadings)-1), 15),
           start.degree=90, clock.wise=FALSE)
# Visualise using a chord diagram!
chordDiagram(bootstrapped_xloadings, 
             # Visualise negative loadings in blue and positive in red
             col = col_fun,
             # Add directionality of the loadings using arrows. 
             directional = 1, direction.type = "arrows", link.arr.type = "big.arrow", diffHeight = -mm_h(2),
             # Scale so that relative contributions are visualised, and only 
             # show the sector names. 
             scale=TRUE, annotationTrack = "name",
             # Sort links according to the width of each sector
             link.sort = TRUE, link.decreasing = TRUE)
# Align vertically. 
abline(v = 0, lty = 2, lwd = 2, col = "#00000080")



### PART 9 - Distribution of Mean Number of Samples Provided by Each AHBA Donor - SUPPLEMENTARY FIGURE 1 ####
# For our main analysis, we used the Schaefer 100-node (17-network) parcellation. 
# This is because whilst we initially wanted to use the Schaefer 400-node parcellation,
# this has a far smaller number of mean samples per AHBA donor for each region,
# thus increasing noise. We justify our decision here.
### PART 9A - Load Data ####
# For each of the 3 parcellations (Schaefer 100-node, Brainnetome 246-node, and 
# Schaefer 400-node, load in their AHBA summary data detailing the number of samples
# each donor provided for each region, and the formatted parcellation labels.)
schaefer100_ahba_counts = read.csv('data/ahba_expression_schaefer100_counts.csv')
# Remove the row index column.
schaefer100_ahba_counts = schaefer100_ahba_counts[,-c(1)]
# Load and format parcellation labels.
schaefer100_labels = read_mat('data/schaefer100x17_1mm_info.mat')
schaefer100_labels = unlist(schaefer100_labels[["schaefer100x17.1mm.info"]][[1]])
schaefer100_labels[1:50] = paste0("lh_",schaefer100_labels[1:50])
schaefer100_labels[51:100] = paste0("rh_",schaefer100_labels[51:100])
# Calculate the number of samples and donors per region.
schaefer100_ahba_counts$'Number of Donors' = rowSums(schaefer100_ahba_counts!=0)
schaefer100_ahba_counts$'Number of Samples' = rowSums(schaefer100_ahba_counts[1:6])
# Add the label column!
schaefer100_ahba_counts$'label' = schaefer100_labels
# Finally, calculate the mean number of samples per donor per region.
schaefer100_ahba_counts$'mean' = schaefer100_ahba_counts$`Number of Samples`/schaefer100_ahba_counts$`Number of Donors`

# And the Brainnetome 246-node parcellation... Note that for this parcellation,
# the MATLAB files provide labels specifying the gyrus of the relevant lobes for
# 246 regions, whilst the ggseg atlas provides the modified cyto-architectonic
# labels for 243 regions. Therefore, we will need to match them up!
brainnetome246_ahba_counts = read.csv('data/ahba_expression_brainnetome246_counts.csv')
brainnetome246_ahba_counts = brainnetome246_ahba_counts[,-c(1)]
brainnetome246_labels = readxl::read_excel("data/BNA_subregions.xlsx")
# For each region, add a 'lh_' or 'rh_' prefix to the cytoarchitectonic label,
# as well as a 'L' or 'R' suffix!
brainnetome246_labels_lh = vector("list",123)
brainnetome246_labels_rh = vector("list",123)
for (ROI in 1:length(brainnetome246_labels_lh)){
  Architectonic_Label = str_split(brainnetome246_labels[ROI,6],",")[[1]][1]
  # Now create a left hemisphere version!
  brainnetome246_labels_lh[[ROI]] = paste0("lh_",Architectonic_Label,"_L")
}
for (ROI in 1:length(brainnetome246_labels_rh)){
  Architectonic_Label = str_split(brainnetome246_labels[ROI,6],",")[[1]][1]
  # Now create a right hemisphere version!
  brainnetome246_labels_rh[[ROI]] = paste0("rh_",Architectonic_Label,"_R")
}
brainnetome246_ahba_counts$'Number of Donors' = rowSums(brainnetome246_ahba_counts!=0)
brainnetome246_ahba_counts$'Number of Samples' = rowSums(brainnetome246_ahba_counts[1:6])
brainnetome246_ahba_counts$'mean' = brainnetome246_ahba_counts$`Number of Samples`/brainnetome246_ahba_counts$`Number of Donors`
brainnetome246_ahba_counts$'label' = c(unlist(brainnetome246_labels_lh), unlist(brainnetome246_labels_rh))

schaefer400_ahba_counts = read.csv('data/ahba_expression_schaefer400_counts.csv')
schaefer400_ahba_counts = schaefer400_ahba_counts[,-c(1)]
schaefer400_labels = read_mat('data/schaefer400x17_1mm_info.mat')
schaefer400_labels = unlist(schaefer400_labels[["schaefer400x17.1mm.info"]][[1]])
schaefer400_labels[1:200] = paste0("lh_",schaefer400_labels[1:100])
schaefer400_labels[201:400] = paste0("rh_",schaefer400_labels[201:400])
schaefer400_ahba_counts$'Number of Donors' = rowSums(schaefer400_ahba_counts!=0)
schaefer400_ahba_counts$'Number of Samples' = rowSums(schaefer400_ahba_counts[1:6])
schaefer400_ahba_counts$'label' = schaefer400_labels
schaefer400_ahba_counts$'mean' = schaefer400_ahba_counts$`Number of Samples`/schaefer400_ahba_counts$`Number of Donors`

### PART 9B - Plot! ####
ggseg(schaefer100_ahba_counts, atlas = schaefer17_100,
      colour = "black", size = 1, mapping = aes(fill = mean),
      position = "stacked") + theme_void() + 
  labs(fill = "Mean Number of \nSamples Per Donor") +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
        legend.key.size = unit(1, 'cm'), legend.key.width = unit(1,'cm'),
        legend.text = element_text(size = 15)) + 
  scale_fill_viridis() +
  guides(fill=guide_colorbar(ticks.colour=NA))

# Note that a warning appears about merhing, but these are usually subcortical 
# regions which can't be visualised on the cortical surface. 
ggseg(brainnetome246_ahba_counts, atlas = brainnetome,
      colour = "black", size = 1, mapping = aes(fill = mean),
      position = "stacked") + theme_void() + 
  labs(fill = "Mean Number of \nSamples Per Donor") +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
        legend.key.size = unit(1, 'cm'), legend.key.width = unit(1,'cm'),
        legend.text = element_text(size = 15)) + 
  scale_fill_viridis() +
  guides(fill=guide_colorbar(ticks.colour=NA))

ggseg(schaefer400_ahba_counts, atlas = schaefer17_400,
      colour = "black", size = 1, mapping = aes(fill = mean),
      position = "stacked") + theme_void() + 
  labs(fill = "Mean Number of \nSamples Per Donor") +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
        legend.key.size = unit(1, 'cm'), legend.key.width = unit(1,'cm'),
        legend.text = element_text(size = 15)) + 
  scale_fill_viridis() +
  guides(fill=guide_colorbar(ticks.colour=NA))

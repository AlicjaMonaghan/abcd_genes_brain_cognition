# This script describes visualization of the distribution of polygenic scores 
# (PGSs) for cognitive ability, for a subset of the (baseline) fourth release 
# of the Adolescent Brain and Cognitive Development (ABCD) study, for which we 
# conducted generative network modelling (GNM). This corresponds to Figure 4. 
# All correspondence to Alicja.Monaghan@mrc-cbu.cam.ac.uk

# Part 1 - Preparing the Work Space 
# Part 2 - Plotting PGS Distribution (Figure 4a)
# Part 3 - Box plot of Eta and Gamma Distributions for Top vs Bottom 10% PGSs (Figure 4b)
# Part 4 - Variability of Nodal Wiring Probabilities (Figure 4c)

### PART 1 - Set Up the Work Space ####
rm(list = ls())
library(ggplot2)
library(viridis)
library(tibble)
library(raveio)
library(dplyr)
# Set to where you saved this directory!
setwd('//cbsu/data/imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_variation_structural_generative_mechanisms_open/data/')
# Load the polygenic scores and GNM parameters.
pgs_and_covariates = read.csv("gnm_parameters_pgs_covariates_anonymised.csv")
# Load the simulated statistics from the 1000 networks generated using the mean 
# eta and gamma for the top and bottom 10% PGS participants, respectively.
low_pgs_simulations = read_mat("group_1_generative_model.mat")
high_pgs_simulations = read_mat("group_2_generative_model.mat")

### PART 2 - Plot PGS Distribution ####
# Create a factor to color the top 10% PGS scores as green, and bottom 10% PGS
# scores as red, with all others as grey. 
nsub = nrow(pgs_and_covariates)
bottom_pgs_boundary = sort(pgs_and_covariates$pgs)[round(nsub*.10)]
top_pgs_boundary = sort(pgs_and_covariates$pgs)[nsub-round(nsub*.10)]
pgs_and_covariates$pgs_level = 
  as.factor(ifelse(pgs_and_covariates$pgs<=bottom_pgs_boundary,"Lowest",
         ifelse(pgs_and_covariates$pgs>=top_pgs_boundary, "Highest","NA")))
# Multiply the PGSs by 10,000 for easier formatting of axis text (avoiding scientific notation).
pgs_and_covariates$pgs = pgs_and_covariates$pgs*10000
# Visualize distribution of polygenic scores
pgs_distribution = ggplot(pgs_and_covariates, mapping = aes(x = pgs, fill = pgs_level)) + 
  geom_histogram(bins=200, alpha = .5) +
  # Ensure to add to the x-axis that the PGSs were multiplied by 10,000 for plotting!
  labs(title = "Distribution of Polygenic Scores for Cognitive Ability",
       x = bquote("Polygenic Scores for Cognitive Ability "(x10^-4)),
       y = "Frequency") +
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15), 
        plot.title = element_text(size=15, face="bold", hjust = .5),
        panel.background = element_blank(), axis.line = element_line(colour="black"),
        axis.ticks = element_blank(), legend.position = "none") +
  # Remove space between histogram bars and x axis
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = c("Lowest" = "#440154FF", "Highest" = "#22A884FF", "NA" = "lightgrey")) 
pgs_distribution
ggsave("pgs_distribution.png", plot = pgs_distribution, device = "png", width=8, height=4, units = "in", dpi = 700)

### PART 3 - Plot GNM Distributions for High vs Low PGS Groups ####
# Create a tibble with the eta and gamma values for participants from the 
# top and bottom 10% tails of the PGS distribution.
low_pgs_indices = sort(pgs_and_covariates$pgs, index.return=T)$ix[1:round(nsub*.10)]
high_pgs_indices = sort(pgs_and_covariates$pgs, index.return=T)$ix[(nsub-round(nsub*.10)+1):nsub]
gnm_parameter_distributions_pgs_contrast = 
  tibble(eta = c(pgs_and_covariates$eta[low_pgs_indices],pgs_and_covariates$eta[high_pgs_indices]),
         gamma = c(pgs_and_covariates$gamma[low_pgs_indices],pgs_and_covariates$gamma[high_pgs_indices]),
         group = factor(rep(c("Lowest","Highest"),each=length(low_pgs_indices))))
# Now create box plots to show the distributions of eta and gamma for each 
# PGS group (top vs bottom 10%). 
gnm_parameter_names = c("eta","gamma")
gnm_parameter_contrast_plot_list = vector("list",2)
for (variable_idx in 1:2){
  variable = gnm_parameter_names[variable_idx]
  gnm_parameter_contrast_plot_list[[variable_idx]] = 
    gnm_parameter_distributions_pgs_contrast %>%
    ggplot(., aes(x = .data[[variable]], y=group, fill=group)) +
    geom_boxplot(colour="grey",alpha=.5, width=.4, outlier.shape = NA, size=1.5) +
    labs(x = parse(text=variable)) +
    theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y = element_blank(),
          axis.ticks.y=element_blank(), axis.text.x = element_text(size=35),
          axis.title.x = element_text(size=60), panel.background = element_blank(),
          axis.line.x = element_line(colour="black"),
          plot.margin = margin(1,1,1,1,"cm"))  +
    scale_fill_manual(values = c("Lowest" = "#440154FF", "Highest" = "#22A884FF")) +
    scale_x_continuous(breaks=round(seq(min(gnm_parameter_distributions_pgs_contrast[[variable]]),
                                  max(gnm_parameter_distributions_pgs_contrast[[variable]]),length.out=5),2))
  print(gnm_parameter_contrast_plot_list[[variable_idx]])
  ggsave(filename = paste0(variable,"_distribution_pgs_contrast.png"),
         device = "png", plot = gnm_parameter_contrast_plot_list[[variable_idx]],
         width=10, height=6, dpi=700)
}

### PART 4 - Visualize Variability of Wiring Probabilities ####
# Following Carozza et al (2023), we shall visualise the frequency of different
# wiring probabilities across time, and the variation in these too! First, 
# average wiring probabilities across all steps, for each group separately.
low_pgs_wiring_probability = low_pgs_simulations[["output/probabilities"]]
high_pgs_wiring_probability = high_pgs_simulations[["output/probabilities"]]
low_pgs_mean_wiring_probability = rowMeans(low_pgs_wiring_probability)
high_pgs_mean_wiring_probability = rowMeans(high_pgs_wiring_probability)
nnetworks = dim(low_pgs_simulations[["output/probabilities"]])[1]
nstep = dim(low_pgs_simulations[["output/probabilities"]])[2]

# For each step, calculate the mean and coefficient of variation (CV) of wiring 
# probabilities.
mean_var_wiring_probabilities_list = vector("list",2)
wiring_probabilities = list(low_pgs_wiring_probability, high_pgs_wiring_probability)
for (group in 1:2){
  # Initialize empty array for this group
  mean_var_wiring_probabilities = array(as.numeric(), dim=c(nnetworks,2))
  # Summary statistics for each network
  for (network in 1:nnetworks){
    mean_var_wiring_probabilities[network,1] = mean(wiring_probabilities[[group]][network,])
    mean_var_wiring_probabilities[network,2] = sd(wiring_probabilities[[group]][network,])/mean(wiring_probabilities[[group]][network,])
  }
  mean_var_wiring_probabilities_list[[group]] = mean_var_wiring_probabilities
}
# Combine both data frames into one!
mean_var_wiring_probabilities_df = 
  bind_rows(data.frame(mean_var_wiring_probabilities_list[[1]]),
            data.frame(mean_var_wiring_probabilities_list[[2]]))
mean_var_wiring_probabilities_df$group = factor(rep(c("Lowest","Highest"),each=nnetworks))
colnames(mean_var_wiring_probabilities_df)[1:2] = c("mean","cv")
# Plot distribution of coefficient of variation!
cv_probabilities = ggplot(mean_var_wiring_probabilities_df,mapping=aes(x=cv,fill=group)) +
  geom_histogram(bins=100,alpha=.5) +
  labs(title = expression(bold("Highest 10% PGS Show Greater Variability \nin Nodal Wiring Probabilities")~bolditalic("P"["i,j"])),
       x = expression("Coefficient of Variation in"~italic("P"["i,j"])), y="Frequency") +
  theme(axis.title = element_text(size=15), axis.text=element_text(size=15),
        plot.title = element_text(size=15, face="bold", hjust=.5),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(), legend.position = "none") +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = c("Lowest" = "#440154FF", "Highest" = "#22A884FF"))
print(cv_probabilities)
ggsave("cv_probabilities.png", plot = cv_probabilities, 
       device = "png", width=8, height=4, units = "in", dpi = 700)

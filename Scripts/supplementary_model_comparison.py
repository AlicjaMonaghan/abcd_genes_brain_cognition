# This script conducts supplementary model comparisons for "Brain wiring economics, network organisation, and population
# level genomics".

from scipy.io import loadmat
import os
import pandas as pd
from scipy.stats import f_oneway, tukey_hsd
import numpy as np
from decimal import Decimal
from bct.utils import adjacency_plot_und

os.chdir('/Users/alicjamonaghan/Desktop/abcd_genes_brain_cognition/')
# Load the group-level generative energy estimates
group_energy_estimates = loadmat('data/group_gnm_energy_99856.mat')['original_simulation_energy']
modeltypes = ['sptl', 'neighbors', 'matching', 'clu-avg', 'deg-avg']
group_energy_estimates_pd = pd.DataFrame(group_energy_estimates, columns=modeltypes)
# Load the seed network used for all simulations and extract the Schaefer 100-node version.
seed = loadmat('data/seed_across_parcellations.mat')['seed']['schaefer100'].item()
# Load the consensus network (group simulation target)
consensus = loadmat('data/group_target_across_parcellations.mat')['target']['schaefer100'].item()
# Find the coordinates for the Schaefer 100-node parcellation
schaefer100_coords = loadmat('data/schaefer100x17_1mm_info.mat', simplify_cells=True)['schaefer100x17_1mm_info']
coords = np.stack([schaefer100_coords['x_mni'], schaefer100_coords['y_mni'], schaefer100_coords['z_mni']], 1)
# Load the generative model parameter estimates and covariates
gnm_parameters_pgs = pd.read_csv('data/gnm_parameters_pgs_covariates.csv')

# Order each column ascendingly
group_energy_estimates_pd = group_energy_estimates_pd.apply(lambda x: x.sort_values().values, axis=0)
# Set the number of simulations across which we're comparing models
n_simulations = [10, 50, 100, 500, 1000]
# Initialise an output array to hold the mean and standard deviation of energy values across wiring rules and number of
# simulations.
model_comparison_descriptive = np.zeros([len(n_simulations), len(modeltypes), 2])
for simulation_index, n_simulation in enumerate(n_simulations):
    # Subset by the number of simulations we're examining. Subtract one to account for zero-indexing.
    energy_subset = group_energy_estimates_pd.loc[0:n_simulation - 1]
    model_comparison_descriptive[simulation_index, :, 0] = energy_subset.mean()
    model_comparison_descriptive[simulation_index, :, 1] = energy_subset.std()
    # Conduct a one-way ANOVA. We expect each ANOVA to be strongly significant.
    anova_statistic, anova_pval = f_oneway(*energy_subset.T.values)
    print('Effect of model type on energy in the top {} simulations is {:.2E}'.format(
        n_simulation, Decimal(anova_pval)))
    # Conduct a post-hoc test (Tukey's honestly significant difference)
    tukey_res = tukey_hsd(*energy_subset.T.values)
    # Extract the comparison between the homophily models

    # print(example.pvalue)
# Save the model comparison descriptive table!
pd.DataFrame(model_comparison_descriptive[:, :, 0], columns=modeltypes).to_csv(
    'data/supplementary_group_model_comparison_mean.csv')
pd.DataFrame(model_comparison_descriptive[:, :, 1], columns=modeltypes).to_csv(
    'data/supplementary_group_model_comparison_std.csv')

# Plot the seed network superimposed onto the consensus
[x, y, z] = adjacency_plot_und((consensus-seed), coords)
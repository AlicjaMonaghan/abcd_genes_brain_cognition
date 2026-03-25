# This script conducts supplementary model comparisons for "Brain wiring economics, network organisation, and population
# level genomics". We also measure the developmental trajectories of global efficiency for simulated networks at each
# step, for those with the top or bottom 10% polygenic scores. Written by Alicja Monaghan in January 2025.

from scipy.io import loadmat
import mat73
import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mode
from bct import edge_nei_overlap_bu, generative_model
from scipy.spatial.distance import pdist


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
# Get the Euclidean distance between coordinates
euclidean = pdist(coords, metric='euclidean')
# Load the generative model parameter estimates, polygenic scores and participant IDs.
optimal_parameters_with_pgs = pd.read_csv('data/ABCD_Individual_Parameters_with_PRS.txt', names=[
    'id', 'eta', 'gamma', 'pc1', 'pgs', 'ethnicity', 'euclidean'], delim_whitespace=True)
# Load the energies associated with the best models and participant IDs.
lowest_energy = mat73.loadmat('data/ABCD_Individual_Lowest_Energy.mat')['ABCD_Individual_Lowest_Energy']
participant_ids = pd.read_csv('data/Processed_Participant_IDs_Accompanying_Sorted_Outputs.txt', header=None)[0]
lowest_energy_pd = pd.DataFrame({'id': participant_ids, 'energy': lowest_energy})
# Load the target connectomes and the seed network for the Schaefer 100-node parcellation
target_connectomes = mat73.loadmat(
    'data/ABCD_Individual_Target_Connectomes_schaefer100x17_Parcellation.mat')['ABCD_Thresholded_27_Streamlines']
seed = mat73.loadmat('data/ABCD_Seed_Network_schaefer100x17_Parcellation.mat')['ABCD_Seed_Network']
# Order each column ascendingly
group_energy_estimates_pd = group_energy_estimates_pd.apply(lambda x: x.sort_values().values, axis=0)
# Set the number of simulations across which we're comparing models
n_simulations = [1, 10, 25, 500, 1000]
# Initialise an output array to hold the mean and standard deviation of energy values across wiring rules and number of
# simulations.
model_comparison_descriptive = np.zeros([len(n_simulations), len(modeltypes), 2])
for simulation_index, n_simulation in enumerate(n_simulations):
    # Subset by the number of simulations we're examining. Subtract one to account for zero-indexing.
    energy_subset = group_energy_estimates_pd.loc[0:n_simulation - 1]
    model_comparison_descriptive[simulation_index, :, 0] = energy_subset.mean()
    model_comparison_descriptive[simulation_index, :, 1] = energy_subset.std()
# Save the model comparison descriptive table!
pd.DataFrame(model_comparison_descriptive[:, :, 0], columns=modeltypes).to_csv(
    'data/supplementary_group_model_comparison_mean.csv')
pd.DataFrame(model_comparison_descriptive[:, :, 1], columns=modeltypes).to_csv(
    'data/supplementary_group_model_comparison_std.csv')

# We observed greater stochasticity in the simulations from participants with the top 10% polygenic scores. We assess
# whether this is a consequence of a weaker fit between the generated and empirical networks for participants in the top
# PGS decile. To test this, we'll compare the minimal energy obtained for participants in each group. First, merge the
# energy estimates with the GNM parameters and PGSs, then sort by increasing PGS.
lowest_energy_pd = pd.merge(lowest_energy_pd, optimal_parameters_with_pgs, on='id', how='right').sort_values(by='pgs')
# Get the number of people in each group
nsub_per_group = int(np.round(len(lowest_energy_pd)*.10))
# Compare the energy using a two-tailed t-test
low_pgs_energy = lowest_energy_pd.head(nsub_per_group).energy
print(f'Energy for low PGS group: mean of {np.round(low_pgs_energy.mean(),3)} and standard deviation of {np.round(low_pgs_energy.std(), 3)}')
high_pgs_energy = lowest_energy_pd.tail(nsub_per_group).energy
print(f'Energy for high PGS group: mean of {np.round(high_pgs_energy.mean(), 3)} and standard deviation of {np.round(high_pgs_energy.std(), 3)}')
energy_comparison_result = ttest_ind(low_pgs_energy, high_pgs_energy)

# Work out how many connections are added from the seed to simulate each child's connectome. Start by removing the
# connectome with index 1145 (note zero-indexing), as it was incorrectly reconstructed.
target_connectomes = np.delete(target_connectomes, 1144, axis=0)
# Create an array which will hold the number of connections to be added for each participant
connections_to_be_added = np.zeros((len(target_connectomes), 1))
for i in range(len(target_connectomes)):
    connections_to_be_added[i] = np.count_nonzero(target_connectomes[i, :, :]) - np.count_nonzero(seed)
print("Adding a minimum of {} connections, maximum of {} and modal {}".format(
    connections_to_be_added.min(), connections_to_be_added.max(), mode(connections_to_be_added)[0][0]))

"""
This script computes two global graph theory metrics - modularity and efficiency - for simulated and empirical
connectomes from the extreme polygenic score deciles from the Adolescent Brain Cognitive Development study. Written by
Alicja Monaghan, MRC Cognition and Brain Sciences Unit, University of Cambridge, in July 2024.
"""

import os
import pandas as pd
import mat73
import numpy as np
from bct import efficiency_bin, modularity_und, assortativity_bin, clustering_coef_bu
from scipy.stats import zscore
import statsmodels.formula.api as sm
from pingouin import mediation_analysis

open_access_dir = '/imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_variation_structural_generative_mechanisms_open'
os.chdir(open_access_dir)
# Import the PGS and optimal GNM parameters
pgs_gnm = pd.read_csv('data/gnm_parameters_pgs_covariates_anonymised.csv', sep=",")
# Add the subject IDs so that we're matching the correct participants...
pgs_gnm_ids = pd.read_csv('data/ABCD_Individual_Parameters_with_PRS.txt', sep=" ",
                          names=['ID', 'eta', 'gamma', 'pc1', 'pgs', 'race', 'euclidean'])
pgs_gnm_complete = pgs_gnm.merge(pgs_gnm_ids, on=['eta', 'gamma', 'pc1', 'pgs'])
# Load the thresholded and binarized structural connectomes
thresholded_structural_connectomes = mat73.loadmat(
    'data/ABCD_Individual_Target_Connectomes_schaefer100x17_Parcellation.mat')['ABCD_Thresholded_27_Streamlines']
schaefer100_participants = pd.read_csv('data/Successfully_Processed_Participants_schaefer100x17_Parcellation.txt')
# Find the indices of participants with structural connectomes and polygenic scores
matching_idx = schaefer100_participants[schaefer100_participants['Participant_IDs'].isin(
    pgs_gnm_complete.ID)].index.tolist()
extracted_structural_connectomes = thresholded_structural_connectomes[matching_idx, :, :]
# Create an array to hold the global graph theory metrics for each thresholded connectome
global_efficiency_array = np.zeros((extracted_structural_connectomes.shape[0]))
modularity_array = np.zeros((extracted_structural_connectomes.shape[0]))
assortativity_array = np.zeros((extracted_structural_connectomes.shape[0]))
global_clustering_array = np.zeros((extracted_structural_connectomes.shape[0]))
nsub = extracted_structural_connectomes.shape[0]
for sub_idx in range(0, nsub):
    global_efficiency_array[sub_idx] = efficiency_bin(extracted_structural_connectomes[sub_idx, :, :])
    modularity_array[sub_idx] = modularity_und(extracted_structural_connectomes[sub_idx, :, :])[1]
    assortativity_array[sub_idx] = assortativity_bin(extracted_structural_connectomes[sub_idx, :, :])
    # Whilst clustering coefficient is a nodal property, we average across nodes to get a global measure
    global_clustering_array[sub_idx] = np.mean(clustering_coef_bu(extracted_structural_connectomes[sub_idx, :, :]))
    print('Processed metrics for {}'.format(pgs_gnm_complete.ID[sub_idx]))
# Create a data frame and save for visualisation!
global_metrics_pd = pd.DataFrame({'global_efficiency': global_efficiency_array, 'clustering': global_clustering_array,
                                  'assortativity': assortativity_array, 'modularity': modularity_array})
global_metrics_pd.to_csv('data/global_graph_theory_metrics.csv', index=False)
# Add the global efficiency values to the PGS table
pgs_gnm_complete['global_efficiency'] = global_efficiency_array
pgs_gnm_complete['modularity'] = modularity_array
lower_bound = np.mean(pgs_gnm_complete.eta) - np.std(pgs_gnm_complete.eta) * 2
upper_bound = np.mean(pgs_gnm_complete.eta) + np.std(pgs_gnm_complete.eta) * 2
# Find participants whose eta values (distance penalty) are more than 2 standard deviations from the mean
filtered_pgs_gnm = pgs_gnm_complete.loc[(pgs_gnm_complete.eta >= lower_bound) & (pgs_gnm_complete.eta <= upper_bound)]
# Standardise the predictors we'll include in our generalised additive mixed models
filtered_pgs_gnm.update(filtered_pgs_gnm[[
    'eta', 'gamma', 'pgs', 'age', 'global_efficiency', 'modularity', 'meanFWD']].apply(zscore))
# Assign rows to deciles based on polygenic scores
filtered_pgs_gnm['pgs_quantile'] = pd.qcut(filtered_pgs_gnm.pgs, 10, labels=False)
# Predict global efficiency and modularity using PGS, mean frame wise displacement, age, and sex.
outcomes = ['global_efficiency', 'modularity']
variables_of_interest = ['eta', 'gamma', 'pgs']
for outcome in outcomes:
    glm_results = sm.glm(formula=outcome + ' ~ pgs + eta + gamma + meanFWD + sex + age', data=filtered_pgs_gnm).fit()
    for variable in variables_of_interest:
        print('Predicting %s using %s: standardised beta of %.3f and p-value of %.3f' % (
            outcome, variable, round(glm_results.params[variable], 3), round(glm_results.pvalues[variable], 3)))
# Predict eta and gamma separately using polygenic scores for cognitive ability, controlling for the same covariates as
# above!
outcomes = ['eta', 'gamma']
for outcome in outcomes:
    glm_results = sm.glm(formula=outcome + ' ~ pgs + meanFWD + sex + age', data=filtered_pgs_gnm).fit()
    print('Predicting %s using polygenic scores: standardised beta of %.3f and p-value of %.3f' %(
        outcome, round(glm_results.params[variable], 3), round(glm_results.pvalues[variable], 3)))
# Now conduct two sets of mediation analyses. Across both, the independent variable is the polygenic scores, and eta is
# the mediator. Global efficiency was the dependent variable in one mediation, whilst modularity was the dependent
# variable in the second analysis.
mediation_analysis(data=filtered_pgs_gnm, x='pgs', m='eta', y='global_efficiency', seed=42, alpha=.05)
mediation_analysis(data=filtered_pgs_gnm, x='pgs', m='gamma', y='global_efficiency', seed=42, alpha=.05)
mediation_analysis(data=filtered_pgs_gnm, x='pgs', m='eta', y='modularity', seed=42, alpha=.05)
mediation_analysis(data=filtered_pgs_gnm, x='pgs', m='gamma', y='modularity', seed=42, alpha=.05)

% This script details finding which genes from the Allen Human Brain Atlas
% (AHBA) significantly predict parameterised nodal wiring costs or values,
% respectively, across 2153 participants from the Adolescent Brain
% Cognitive Development (ABCD) study. For each participant, we ran two
% partial least squares regressions (PLSs), both of which used AHBA genes
% (1 x 12431 matrix) as predictors. The outcome in the first PLS is
% parameterised nodal wiring costs, and the outcome in the second PLS is
% parameterised nodal wiring values. You can run these PLSs for each
% participant using the run_pls_regression_ahba_nodal_parameters.mat code.
% Else, you can load our outputs here. Note that for the PLSs, we used
% 10,000 bootstraps to determine confidence intervals, and 10,000
% permutations to create null distributions. Correspondence to Alicja
% Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk

% STEPS:
% 1. Set up work space and load data.
% 2. Find AHBA genes predictive of parameterised nodal wiring costs, with a
% permuted p-value of less than .05, and sort by descending X loadings onto
% the first PLS component.
% 3. Repeat step 2 but using parameterised nodal wiring value.

%% PART 1 - Set up the workspace.
clear;clc;
% SET WORKING DIRECTORY TO WHERE YOUR SAVED DIRECTORY IS!
cd('abcd_genomic_variation_structural_generative_mechanisms_open/');
% Load the AHBA gene expression arrays.
ahba_gene_expression = readtable('data/AHBA_expression_rnaseq_schaefer100_cleaned.csv');
% Retrieve the probe names associated with each gene
ahba_gene_names = ahba_gene_expression.Properties.VariableNames';
% Remove the label column.
ahba_gene_names(1,:) = [];
% And convert to a string.
ahba_gene_names = string(ahba_gene_names);

% Load the permuted probability coefficient for each gene, X loadings, and 
% variance explained in the PLSs with parameterised nodal wiring costs as 
% the outcome. nodal_costs_pcorr is an nsub x ncomp x ngenes array.
nodal_costs_pcorr = load('data/parameterised_nodal_costs_pls_permuted_pcorr.mat');
nodal_costs_pcorr = nodal_costs_pcorr.Nodal_Costs_PLS_Pcorr;                                      
% nodal_costs_xloading is an nsub x ngenes x ncomp array.
nodal_costs_xloading = load('data/parameterised_nodal_costs_pls_xloadings.mat');
nodal_costs_xloading = nodal_costs_xloading.Nodal_Costs_PLS_Xloading;
% nodal_costs_var is an nsub x ncomp x 2 array.
nodal_costs_var = load('data/parameterised_nodal_costs_pls_variance.mat');
nodal_costs_var = nodal_costs_var.Nodal_Costs_PLS_var;

% Do the same for the PLSs with parameterised nodal wiring value as the
% outcome!
nodal_value_pcorr = load('data/parameterised_nodal_value_pls_permuted_pcorr.mat');
nodal_value_pcorr = nodal_value_pcorr.Nodal_Value_PLS_Pcorr;
nodal_value_xloading = load('data/parameterised_nodal_value_pls_xloadings.mat');
nodal_value_xloading = nodal_value_xloading.Nodal_Value_PLS_Xloading;
nodal_value_var = load('data/parameterised_nodal_value_pls_variance.mat');
nodal_value_var = nodal_value_var.Nodal_Value_PLS_var;

%% PART 2 - Finding AHBA Genes Predictive of Parameterised Nodal Wiring Costs 
% We shall find which AHBA genes are significant predictors of 
% parameterised nodal costs, rank them according to their loading onto the 
% first principal component.

% Find the % variance of nodal costs explained by AHBA across participants, 
% for each component. We find that the SECOND component accounts for the
% majority of variance (~38.6%). However, we shall extract genes using the
% first principal component, as the genetic data set is larger than the 
% parameters data set. When we plot the loadings, we see that the first
% component captures an anterior-posterior axis. 
mean_variance_explained_nodal_costs_onto_ahba = mean(squeeze(nodal_costs_var(:,:,2)));
sd_variance_explained_nodal_costs_onto_ahba = std(squeeze(nodal_costs_var(:,:,2)));

% Now we shall find which genes were significant predictors of
% parameterised nodal wiring costs across which proportion of participants. 
% Again, we shall use the first principal component as our reference. For 
% each participant, loop through the p-value vector and keep the indices of
% the genes with a permuted p-value of less than .05
nsub = size(nodal_costs_pcorr,1);
nodal_costs_predictive_genes_idx = cell(nsub, 1);
for sub = 1:nsub
    [row,~] = find(squeeze(nodal_costs_pcorr(sub,1,:)) < .05);
    % Find the row indexing genes predictive of parameterised nodal wiring
    % costs for this participant, with a permuted p value of less than .05.
    nodal_costs_predictive_genes_idx{sub,1} = row;
end
% Find the unique indices/genes.
nodal_costs_predictive_genes_unique = unique(cell2mat(nodal_costs_predictive_genes_idx),'rows');
% Average PLS loadings for the first component across all participants. 
mean_xloading_nodal_costs_first_component = mean(nodal_costs_xloading(:,:,1),1);
% And find the average PLS loadings for the genes which predict
% parameterised nodal wiring costs, and convert to a table.
nodal_costs_predictive_genes_unique_table = table(mean_xloading_nodal_costs_first_component(nodal_costs_predictive_genes_unique)');
% Format the table.
nodal_costs_predictive_genes_unique_table.Properties.VariableNames = {'PLS1_Loading'};
% Extract the associated AHBA gene names.
nodal_costs_predictive_genes_unique_table.('ahba_gene_names') = ahba_gene_names(nodal_costs_predictive_genes_unique);
% Sort by decreasing gene loading onto PLS1
nodal_costs_predictive_genes_unique_table = sortrows(nodal_costs_predictive_genes_unique_table,1,'descend');
% And update user...
fprintf('%d genes predict parameterised nodal wiring costs.\n',height(nodal_costs_predictive_genes_unique_table));
writetable(nodal_costs_predictive_genes_unique_table,'data/parameterised_nodal_wiring_costs_predictive_genes.txt','Delimiter',' ');
%% PART 3 - Repeat Sorting Process for Parameterised Nodal Wiring Value!
% Find the % variance of nodal costs explained by AHBA across participants, 
% for each component. We find that the FIRST component accounts for the
% majority of variance (25.9%). 
mean_variance_explained_nodal_value_onto_ahba = mean(squeeze(nodal_value_var(:,:,2)));
sd_variance_explained_nodal_value_onto_ahba = std(squeeze(nodal_value_var(:,:,2)));
nodal_value_predictive_genes_idx = cell(nsub, 1);
for sub = 1:nsub
    [row,~] = find(squeeze(nodal_value_pcorr(sub,1,:)) < .05);
    % Extract the associated AHBA names and assign to output matrix.
    nodal_value_predictive_genes_idx{sub,1} = row;
end
% Find the unique indices/genes.
nodal_value_predictive_genes_unique = unique(cell2mat(nodal_value_predictive_genes_idx),'rows');
% Average PLS loadings for the first component across all participants. 
mean_xloading_nodal_value_first_component = mean(nodal_value_xloading(:,:,1),1);
% And find the average PLS loadings for the genes which predict
% parameterised nodal wiring value, and convert to a table.
nodal_value_predictive_genes_unique_table = table(mean_xloading_nodal_value_first_component(nodal_value_predictive_genes_unique)');
% Format the table.
nodal_value_predictive_genes_unique_table.Properties.VariableNames = {'PLS1_Loading'};
% Extract the associated AHBA gene names.
nodal_value_predictive_genes_unique_table.('ahba_gene_names') = ahba_gene_names(nodal_value_predictive_genes_unique);
% Sort by decreasing gene loading onto PLS1
nodal_value_predictive_genes_unique_table = sortrows(nodal_value_predictive_genes_unique_table,1,'descend');
fprintf('%f genes predict parameterised nodal wiring value.\n',height(nodal_value_predictive_genes_unique_table));
% And save...
writetable(nodal_value_predictive_genes_unique_table,'data/parameterised_nodal_wiring_value_predictive_genes.txt','Delimiter',' ');
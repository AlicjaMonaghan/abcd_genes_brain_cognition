% This script details simulating matching-neighbours generative network
% models (GNMs) for mean eta and gamma derived from participants with the
% bottom and top 10% polygenic scores for cognitive ability. We test
% for significant differences in graph theory metrics, such as efficiency
% and modularity. 

%% PART 1 - Set up the workspace.
% Load up the ABCD individual eta and gamma parameters with polygenic
% scores, and format! Convert the variables to numeric as appropriate.
clear;clc;
cd('/Imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_variation_structural_generative_mechanisms_open/');
% Load GNM parameters (eta and gamma), polygenic scores for general
% cognitive ability, cognitive ability scores, and covariates.
parameters_and_covariates = readtable('data/gnm_parameters_pgs_covariates_anonymised.csv');
% Extract the individual variables from the data.
eta = parameters_and_covariates.eta;
gamma = parameters_and_covariates.gamma;
pgs = parameters_and_covariates.pgs;
% TO WHERE YOU HAVE SAVED THIS DIRECTORY
% Add path folder with toolboxes and functions. This script requires the 
% brain connectivity toolbox (Rubinov and Sporns, 2010) and a function to
% evaluate GNM fit through energy (fcn_ks).
addpath('/imaging/astle/users/da04/PhD/toolboxes/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
% Load the coordinates for the Schaefer 100-node parcellation, and
% calculate the Euclidean distance between points. 
schaefer100x17_1mm_info = load('data/schaefer100x17_1mm_info.mat');
schaefer100x17_1mm_info = schaefer100x17_1mm_info.schaefer100x17_1mm_info;
parcellation_coordinates = [schaefer100x17_1mm_info.x_mni,...
            schaefer100x17_1mm_info.y_mni, schaefer100x17_1mm_info.z_mni];
D = squareform(pdist(parcellation_coordinates));
% For the GNMs, set the seed, number of regions of interest, model
% variables, 
seed_struct = load('data/seed_across_parcellations.mat');
seed = seed_struct.seed.schaefer100;
nroi = size(seed,1);
modelvar = [{'powerlaw'},{'powerlaw'}]; 
% For the GNM, we need to simulate X number of connections. To allow valid
% comparisons between PGS groups, we need to simulate the same number of
% connections. Thus, use the group-level connectome as target.
target_struct = load('data/group_target_across_parcellations.mat');
target = target_struct.target.schaefer100;
% Set number of bi-directional connections
m = nnz(target)/2;
% And network cardinality
n = length(target);

%% PART 2 - Testing for Differences in Eta and Gamma by PGS Groups %%
% Sort the polygenic scores decreasingly, and select the top and bottom 10%
[sorted_pgs,sorted_pgs_idx]= sort(pgs);
% The first and second slices of the third dimension holds the bottom and
% top 10% percentiles of polygenic scores, respectively.
tail_end = .1;
variables = [eta, gamma, pgs];
percentile_array = zeros([round(length(eta)*tail_end),length(variables),2]);
% Set the number of participants in each group
nsub = round(length(eta)*tail_end);
for variable_idx = 1:3
    % Assign the bottom 10% percentile for each measure
    percentile_array(:,variable_idx,1) = variables(sorted_pgs_idx(1:nsub),variable_idx);
    % And the top 10% percentile
    percentile_array(:,variable_idx,2) = variables(sorted_pgs_idx(end-nsub:end-1),variable_idx);
end
% Test for a significant difference between the percentiles in terms of eta
% and gamma. Calculate the mean parameters per group (bottom vs top 10%). 
gnm_parameter_names = ["eta","gamma"];
for variable_idx = 1:2
    [h,p,ci,stats] = ttest2(squeeze(percentile_array(:,variable_idx,1)),squeeze(percentile_array(:,variable_idx,2)));
    effect_size = meanEffectSize(squeeze(percentile_array(:,variable_idx,1)), squeeze(percentile_array(:,variable_idx,2)));
    if p < .05
        fprintf(['Two-sample t-test revealed a significant difference between the ' ...
            'means of the top and bottom %.d percent PGSs in terms of %s, ' ...
            'where t(%d) = %.3f, p = %.3f\n'], ...
            tail_end*100,gnm_parameter_names(variable_idx),stats.df,stats.tstat,p);
    else
        fprintf(['Two-sample t-test revealed no significant difference between the ' ...
            'means of the top and bottom %.d percent PGSs in terms of %s, ' ...
            'where t(%d) = %.3f, p = %.3f\n'], ...
            tail_end*100,gnm_parameter_names(variable_idx),stats.df,stats.tstat,p);
    end
    fprintf(['Effect size is %.3f, with 95 percent confidence intervals of %.3f and %.3f.\n'], ...
        round(effect_size.Effect, 3), round(effect_size.ConfidenceIntervals(1), 3), ...
        round(effect_size.ConfidenceIntervals(2), 3));
end

%% PART 3 - Simulate Structural Connectivity for Bottom and Top 10% PGS Groups
pgs_groups = ["bottom 10%", "top 10%"];
% All graph theory measures collected are NODAL.
nodal_graph_theory_measures = ["local efficiency", "clustering", "participation coefficient"];
global_graph_theory_measures = ["global efficiency", "modularity"];
nsim = 1000;
% Create an output array for the graph theory metrics for the simulations
% of both PGS groups!
y_nodal = zeros([length(pgs_groups), length(nodal_graph_theory_measures), nroi, nsim]);
y_global = zeros([length(pgs_groups), length(global_graph_theory_measures), nsim]);
% And an array for the mean nodal graph theory metrics across simulations
y_nodal_mean = zeros([length(pgs_groups),length(nodal_graph_theory_measures),nroi]);
for pgs_group_idx = 1:length(pgs_groups)
    % Mean eta and gamma for PGS group...
    mean_eta = mean(squeeze(percentile_array(:,1,pgs_group_idx)));
    mean_gamma = mean(squeeze(percentile_array(:,2,pgs_group_idx)));
    fprintf('The %s PGS group has mean eta of %.3f and mean gamma of %.3f\n', ...
        pgs_groups(pgs_group_idx),mean_eta,mean_gamma);
    % Create 1000 simulated matching-neighbours networks using these means.
    for sim = 1:nsim
        B = generative_model(seed,D,m,"neighbors",modelvar,[mean_eta, mean_gamma]);
        nB = size(B,2);
        b = zeros(n); 
        b(B(:,1)) = 1; 
        b = b + b'; 
        % Calculate each graph theory measure, starting with nodal.
        y_nodal(pgs_group_idx,1,:,sim) = efficiency_bin(b);
        y_nodal(pgs_group_idx,2,:,sim) = clustering_coef_bu(b);
        y_nodal(pgs_group_idx,3,:,sim) = participation_coef(b,0);
        % And now global, starting with global efficiency.
        y_global(pgs_group_idx,1,sim) = efficiency_bin(b);
        [~, Q] = modularity_und(b);
        y_global(pgs_group_idx,2,sim) = Q;
        clear B b Q
    end
    % After collecting graph theory metrics for all simulations for this
    % PGS group, average across simulations.
    for measure = 1:length(nodal_graph_theory_measures)
        for node = 1:nroi
            y_nodal_mean(pgs_group_idx,measure,node) = mean(squeeze(y_nodal(pgs_group_idx,measure,node,:)));
        end
    end
    clear mean_eta mean_gamma
end
% Conduct a two-sample t-test comparing PGS groups in terms of graph theory
% metrics for their simulated networks.
for measure = 1:length(nodal_graph_theory_measures)
    [~,p,~,stats] = ttest2(squeeze(y_nodal_mean(1,measure,:)), squeeze(y_nodal_mean(2,measure,:)));
    % Also calculate Cohen's d as an effect size 
    effect_size = meanEffectSize(squeeze(y_nodal_mean(1,measure,:)),squeeze(y_nodal_mean(2,measure,:)));
    if p < .05
        fprintf(['Significant difference between the means of the top and ' ...
            'bottom 10 percent PGSs in terms of %s, where t(%d) = %.3f, p = %.3f\n'], ...
            nodal_graph_theory_measures(measure),stats.df,stats.tstat,p);
    else
        fprintf(['No significant difference between the means of the top and ' ...
            'bottom 10 percent PGSs in terms of %s, where t(%d) = %.3f, p = %.3f\n'], ...
            nodal_graph_theory_measures(measure),stats.df,stats.tstat,p);
    end
    fprintf('Effect size for %s is %.4f, with CIs of %.4f and %.4f, respectively.\n', ...
        nodal_graph_theory_measures(measure),effect_size.Effect, ...
        effect_size.ConfidenceIntervals(1), ...
        effect_size.ConfidenceIntervals(2));
    clear p stats
end
% Now test for significant differences in global graph theory metrics, such
% as global efficiency and modularity.
for measure = 1:length(global_graph_theory_measures)
    [~,p,~,stats] = ttest2(squeeze(y_global(1,measure,:)), squeeze(y_global(2,measure,:)));
    effect_size = meanEffectSize(squeeze(y_global(1,measure,:)), squeeze(y_global(2,measure,:)));
    if p < .05
        fprintf(['Significant difference between the means of the top and ' ...
            'bottom 10 percent PGSs in terms of %s, where t(%d) = %.3f, p = %.3f\n'], ...
            global_graph_theory_measures(measure),stats.df,stats.tstat,p);
    else
        fprintf(['No significant difference between the means of the top and ' ...
            'bottom 10 percent PGSs in terms of %s, where t(%d) = %.3f, p = %.3f\n'], ...
            global_graph_theory_measures(measure),stats.df,stats.tstat,p);
    end
    fprintf('Effect size for %s is %.4f, with CIs of %.4f and %.4f, respectively.\n', ...
        global_graph_theory_measures(measure),effect_size.Effect, ...
        effect_size.ConfidenceIntervals(1), ...
        effect_size.ConfidenceIntervals(2));
    clear p stats
end



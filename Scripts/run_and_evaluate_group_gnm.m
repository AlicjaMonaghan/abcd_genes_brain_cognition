% This script details running and evaluating the performance of 5 wiring
% rules in group-level generative network models (GNMs). See Vertes et al., 
% 2012, and Betzel et al., 2016 for further details about GNMs. These rules
% represent the best-performing (based on prior work) out of 4 GNM classes:
% spatial ("sptl", edges added between spatially proximal nodes), neighbors
% (a homohpily rule where edges are added between nodes with similar 
% neighbours), matching (a homophily rule where edges are added between 
% nodes with similar connectivity profiles), average clustering ("clu-avg",
% edges added between nodes with high average clustering coefficients), and
% average degree ("deg-avg", edges added between highly connected nodes). 
% Note we evaluated 2 rules from the homophily class, as this class has
% consistently performed the best in previous work. 

% STEPS:
% STEP 1: Preparing the work space
% STEP 2: Initialising GNMs across 3 parcellations: Schaefer 100-node,
% Brainnetome 246-node, and Schaefer 400-node. 
% STEP 3: Retrieving statistics from observed connectomes. 
% STEP 4: Running GNMs across 99,856 unique eta-gamma combinations. Note 
% that you can either run the GNMs yourself, or load the data in the code 
% with the original energy values from our simulations (which we provide 
% for the Schaefer 100-node parcellation). Also, even if you use the same 
% number of runs and parameter combinations as our study, your results may 
% vary slightly due to GNMs being probabilistic. 
% STEP 5: Extracting simulated statistics from 1000 simulations of each
% model's lowest-energy parameter combination.
% STEP 6: Comparing GNM rules through topological dissimilarity, across
% parcellations and the 1000 lowest-energy parameter simulations. For each
% parcellation, we conduct a one-way ANOVA between models, and post-hoc
% tests. This forms SUPPLEMENTARY TABLE 7. See Akarca et al., 2022
% (https://www.biorxiv.org/content/10.1101/2022.03.09.483605v1.full) for
% further details.
% STEP 7: Repeat step 6, but with correlations of degree from simulated
% connectomes with observed connectomes. 
% STEP 8: Visualising energy landscapes - FIGURE 2D
% STEP 9: Comparing energy from top 25 lowest-energy simulations for the
% Schaefer 100-node parcellation.
% STEP 10: Visualising topological dissimilarity fingerprints - FIGURE 2E.
% Note that any run can be visualised in order to represent the general
% statistical trends. 

% This code was based off of work by Dr. Danyal Akarca, at the MRC
% Cognition and Brain Sciences Unit, University of Cambridge.
% Correspondence to Alicja Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk

%% PART 1 - PREPARE THE WORKSPACE %%
clear;clc;
cd(['//cbsu/data/imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_' ...
    'variation_structural_generative_mechanisms_open/']); % <--- SET THIS 
% TO WHERE YOU HAVE SAVED THIS DIRECTORY
% Add path folder with toolboxes and functions. This script requires the 
% brain connectivity toolbox (Rubinov and Sporns, 2010) and a function to
% evaluate GNM fit through energy (fcn_ks).
addpath('toolboxes_and_functions/');
addpath('toolboxes_and_functions/2019_03_03_BCT/');
% We will loop over 3 parcellations: Schaefer 100-node (17-network),
% Brainnetome 246-node, and Schaefer 400-node (17-network). For each, we
% shall extract the lowest-energy parameter combination, conduct 1000
% simulations for each rule's lowest-energy parameter combination, test for
% significant differences in topological dissimilarity and capturing
% observed degree between model classes!
parcellation_list = ["schaefer-100", "brainnetome-246","schaefer-400"];
% Load the seed and target structures.
seed_struct = load('data/seed_across_parcellations.mat');
target_struct = load('data/group_target_across_parcellations.mat');
% Create a cell to hold graph theory metrics of the observed connectomes
total_number_of_comparison_measures = 7;
x = cell(length(parcellation_list),total_number_of_comparison_measures); 
% Initialise the GNM rules. If you want to run your own simulations across
% the original 99856 parameter combinations, initialise nparams too.
modeltypes = ["sptl", "neighbors", "matching", "clu-avg", "deg-avg"];
nparams = 99865;
% Initialise output arrays for the GNMs.
group_structural_gnm_output = struct;
group_structural_gnm_output.ks = zeros(length(parcellation_list),length(modeltypes),nparams,4);
% Set the number of repeats we'll use for the lowest-energy parameter
% combination for each rule and parcellation.
nrep = 1000;
% Load the lowest-energy eta and gamma combinations for each parcellation
% (rows) and model (columns).
lowest_energy_eta = load('data/lowest_energy_eta_group.mat');
lowest_energy_eta = lowest_energy_eta.Best_Eta;
lowest_energy_gamma = load('data/lowest_energy_gamma_group.mat');
lowest_energy_gamma = lowest_energy_gamma.Best_Gamma;

%% PART 2 - INITIALISE GNMs %%
% Here, we loop over each of 3 parcellations and retrieve observed 
% statistics. In Part 4, you can run your own GNMs using the full 99865
% parameter combinations. In Part 5, we conduct 1000 simulations of each
% model's lowest-energy parameter combination, and calculate topological
% dissimilarity and correlation with observed degree, which is tabulated in
% Supplementary Table 7. 
for parcellation = 1:3
    % Update user...
    fprintf('Processing the %s parcellation.\n',parcellation_list(parcellation));
    % Load the seed and target structures.
    if parcellation == 1
        % Load the parcellation meta data. 
        parcellation_metadata = load('data/schaefer100x17_1mm_info.mat');
        % Load the coordinates for the parcellation
        parcellation_coordinates = [parcellation_metadata.schaefer100x17_1mm_info.x_mni,...
            parcellation_metadata.schaefer100x17_1mm_info.y_mni,...
            parcellation_metadata.schaefer100x17_1mm_info.z_mni];
        % Set the seed.
        seed = seed_struct.seed.schaefer100;
        % Set the GNM target.
        target = target_struct.target.schaefer100;
    elseif parcellation == 2
        parcellation_metadata = load('data/brainnetome246_info.mat');
        % Load the coordinates for the parcellation
        parcellation_coordinates = [parcellation_metadata.brainnetome246.x_mni,...
            parcellation_metadata.brainnetome246.y_mni,...
            parcellation_metadata.brainnetome246.z_mni];
        % Set the seed.
        seed = seed_struct.seed.brainnetome246;
        % Set the GNM target.
        target = target_struct.target.brainnetome246;
    else
        parcellation_metadata = load('data/schaefer400x17_1mm_info.mat');
        % Load the coordinates for the parcellation
        parcellation_coordinates = [parcellation_metadata.schaefer400x17_1mm_info.x_mni,...
            parcellation_metadata.schaefer400x17_1mm_info.y_mni,...
            parcellation_metadata.schaefer400x17_1mm_info.z_mni];
        % Set the seed
        seed = seed_struct.seed.schaefer400;
        % Set the GNM target.
        target = target_struct.target.schaefer400;
    end
    % Set up an output array for the simulated statistics for this
    % parcellation.
    nroi = size(seed,2);
    y = zeros(length(modeltypes),total_number_of_comparison_measures,nroi,nrep);
    % Calculate the euclidean distance between the coordinates.
    D = squareform(pdist(parcellation_coordinates));
    % Set the number of regions in the parcellation
    nroi = size(D,1);
    % These are the 5 wiring rules to be assessed.
    modeltypes = ["sptl", "neighbors", "matching", "clu-avg", "deg-avg"];
    modelvar = [{'powerlaw'},{'powerlaw'}]; 
    % Set up eta and gamma limits. These were selected to be as wide as
    % possible to ensure we find the lowest-energy combination.
    etalimits = [-7,7];
    gamlimits = [-7,7];
    % Set the number of runs for the model i.e. eta-gamma combinations. Note
    % that in our study, we used 100,000 runs.
    nruns = 100000;
    % Here we are initialising a grid search with eta-gamma combinations 
    % bounded by the limits above. spanning the eta and gamma limits,
    [p,q] = meshgrid(linspace(etalimits(1),etalimits(2),sqrt(nruns)), ...
        linspace(gamlimits(1),gamlimits(2),sqrt(nruns)));
    params = unique([p(:) q(:)],'rows');
    nparams = size(params,1);
    %% PART 3 - Retrieve Statistics from Observed Connectomes %%
    % Get key observed statistics. Note that "KS" denotes Kolmogorov-Smirnov
    % test. The model with the lowest energy is the best fit, where energy is
    % the largest of 4 KS statistics comparing the distributions of nodal 
    % degree, total edge length, betweenness-centrality and clustering 
    % coefficients between empirical connectomes and those simulated by the
    % GNM. 

    % The 'target' connectome is what we are trying to model!
    % Number of bi-directional connections
    m = nnz(target)/2; 
    % Network cardinality.
    n = length(target); 
    % Find key graph theory properties for the empirical (observed) connectome.
    % Sum the number of connections in each column (number of streamlines), 
    % hence DEGREE
    x{parcellation,1} = sum(target,2); 
    % Clustering coefficient
    x{parcellation,2} = clustering_coef_bu(target); 
    % Betweenness-centrality
    x{parcellation,3} = betweenness_bin(target)'; 
    % Find the upper triangle of the connectome larger than 0, and make into a 
    % double. This is edge length.
    x{parcellation,4} = D(triu(target,1) > 0); 
    % The following are additional metrics not included in the energy equation
    % (which tells us the best-fitting model), but will be used to test the
    % ability of the GNMs to recapitulate properties not included in the energy
    % equation, and other measures of fit, such as topological dissimilarity. 
    % Local efficiency.
    x{parcellation,5} = efficiency_bin(target,1);
    % Eigenvector centrality. 
    x{parcellation,6} = eigenvector_centrality_und(target);
    % Modularity - Newman's (2006) spectral community detection algorithm
    x{parcellation,7} = modularity_und(target);
    
    %% PART 4 - Option to Run Your Own GNMs Across 99,856 Parameter Combinations %%
    % If you want to run your own simulations, do so here: 
    % group_structural_gnm_output.energy = zeros(length(modeltypes),nparams);
    % for model = 1:length(modeltypes)
    %     y = zeros(length(modeltypes),total_number_of_comparison_measures,nroi);
    %     fprintf('Running the %s model.\n',modeltypes{model});
    %     B = generative_model(seed,D,m,modeltypes{model},modelvar,params);
    %     nB = size(B,2);
    %     % Create a matrix of zeros for each measure of model fit
    %     K = zeros(nB,comparison_measures_for_ks); 
    %     % For each eta-gamma combination, calculate model energy.
    %     for iB = 1:nB
    %         % Create a matrix of zeros for the cardinality of the network.
    %         b = zeros(n); 
    %         % For iB indexing the parameter of the generative network model, 
    %         % find the associated column in the cardinality network. 
    %         b(B(:,iB)) = 1;
    %         % Represent b as an adjacency matrix 
    %         b = b + b'; 
    %         % Sum the number of connections in each column (number of
    %         % streamlines). This is degree.
    %         y(model,1,:) = sum(b,2); 
    %         % Clustering coefficient.
    %         y(model,2,:) = clustering_coef_bu(b); 
    %         % Betweenness-centrality.
    %         y(model,3,:) = betweenness_bin(b)'; 
    %         % Edge length.
    %         y(model,4,:) = D(triu(b,1) > 0); 
    %         % Calculate the energy of the synthetic network using degree, 
    %         % clustering coefficient, betweenness centrality, and edge length 
    %         % distributions (see Betzel et al., 2016).
    %          for j = 1:comparison_measures_for_ks 
    %              %KS function as defined by Betzel et al., 2016
    %              K(iB,j) = fcn_ks(x{parcellation,model,j},y(model,i,:)); 
    %          end
    %         % Keep the energy for this model and eta-gamma specification. 
    %         group_structural_gnm_output.energy(model,iB) = max(K(iB,:));
    %     end
    % end
    % 

    %% PART 5 - EXTRACT SIMULATED LOCAL STATISTICS FROM 1000 REPEATS OF LOWEST-ENERGY GNM PARAMETER COMBINATIONS %%
    % For each model, find the lowest energy simulation, and extract the
    % synthetic connectome associated with the lowest-energy eta-gamma
    % combination. We repeat this 1000 times to find an overall trend of which 
    % model performs best. This is because GNMs are probabilistic, and so the 
    % results will vary slightly for each simulation of the same parameter 
    % combination.
    
    % Repeat the GNM 1000 times for each model's lowest-energy parameter
    % combination!
    for t = 1:nrep
        % Loop through each of the models, and re-run the simulations for
        % the lowest-energy eta and gamma combination. 
        for model = 1:length(modeltypes)
            optimal_params = [lowest_energy_eta(parcellation,model), lowest_energy_gamma(parcellation,model)];
            % Re-run the simulation for the optimal parameters. 
            B = generative_model(seed,D,m,modeltypes{model},modelvar,optimal_params);
            % Calculate graph theory metrics for this simulation. We'll compare
            % these with the metrics of the observed connectome, through
            % topological dissimilarity and spatial layout of nodes, to really
            % ensure that we've selected the best model. Assign to the output.
            nB = size(B,2);
            for iB = 1:nB
                b = zeros(n);
                b(B(:,iB)) = 1;
                b = b + b';
                % Degree. 
                y(model,1,:,t) = sum(b,2); 
                % Clustering coefficient.
                y(model,2,:,t) = clustering_coef_bu(b); 
                % Betweenness-centrality.
                y(model,3,:,t) = betweenness_bin(b)'; 
                % Find the total edge length, for each node! Note that this differs
                % from the energy calculation, which would just use edge length.
                edge_length = D(triu(b,1) > 0); 
                total_edge_length = zeros(nroi,1);
                for n = 1:nroi
                    s = edge_length(n,:);
                    y(model,4,n,t) = sum(norm(n-s));
                end
                % Collect the 3 additional nodal metrics not included in the energy
                % equation: local efficiency, eigenvector centrality, and
                % modularity.
                y(model,5,:,t) = efficiency_bin(b,1);
                y(model,6,:,t) = eigenvector_centrality_und(b);
                y(model,7,:,t) = modularity_und(b);
            end
                fprintf('Extracted simulated statistics for the group-level %s model for run %f.\n',modeltypes(model),t);
        end
    end
    % Save the simulated statistics for this parcellation!
    filename = sprintf('data/simulated_statistics_lowest_energy_eta_gamma_1000_runs_%s.mat',parcellation_list(parcellation));
    save(filename,"y");
    fprintf('Saved simulated statistics for the %s parcellation!\n',parcellation_list(parcellation));
    clear y
end
%% PART 6 - COMPARE GNM RULES THROUGH TOPOLOGICAL DISSIMILARITY %%
% Topological dissimilarity tests to what extent the simulated connectome
% captures the LOCAL distribution of graph theory metrics as the observed
% connectome. Note that these graph theory metrics are not only those
% included in the energy equation, hence testing generalisability of
% performance. The smaller the topological dissimilarity, the better the
% model fit! 

% Load the observed statistics for each of the parcellations.
x = load('data/observed_statistics_group_gnm_all_parcellations.mat');
x = x.x;

for parcellation = 1:length(parcellation_list)
    if parcellation == 1
        nroi = 100;
        % Load the simulated statistics from 1000 runs of the lowest-energy
        % parameter combination for each model for this parcellation.
        simulated_stats = load('data/simulated_statistics_lowest_energy_eta_gamma_1000_runs_schaefer-100.mat');
        simulated_stats = simulated_stats.y;
    elseif parcellation == 2
        nroi = 246;
        simulated_stats = load('data/simulated_statistics_lowest_energy_eta_gamma_1000_runs_brainnetome-246.mat');
        simulated_stats = simulated_stats.y;
    else 
        nroi = 400;
        simulated_stats = load('data/simulated_statistics_lowest_energy_eta_gamma_1000_runs_schaefer-400.mat');
        simulated_stats = simulated_stats.y;
    end

    % Initialise output arrays for topological dissimilarity and
    % correlations with observed nodal degree across 1000 simulations of
    % the lowest-energy parameter combinations for each model.
    dissimilarity_repeat = zeros(length(modeltypes),nrep);
    spatial_repeat = zeros(length(modeltypes),nrep);
    % Create an array which will hold all statistics, both observed and
    % simulated, for this parcellation.
    statistics = zeros(1+length(modeltypes),nroi,total_number_of_comparison_measures);
    % Assign the observed statistics.
    statistics(1,:,1) = x{parcellation,1};
    statistics(1,:,2) = x{parcellation,2};
    statistics(1,:,3) = x{parcellation,3};
    statistics(1,:,5) = x{parcellation,5};
    statistics(1,:,6) = x{parcellation,6};
    statistics(1,:,7) = x{parcellation,7};

    % For topological dissimilarity, calculate the total distance from each
    % node.
    total_nodal_edge_length_empirical = zeros(nroi,1);
    for n = 1:nroi
        s = x{parcellation,4}(n,:);
        total_nodal_edge_length_empirical(n) = sum(norm(n-s));
    end
    statistics(1,:,4) = total_nodal_edge_length_empirical';
        
    % Initialise an output variable for correlations between simulated and
    % observed nodal metrics for each model.
    correlated_simulated_observed_across_models = zeros(1+length(modeltypes), ...
        total_number_of_comparison_measures,total_number_of_comparison_measures);

    for t = 1:nrep
        % Assign the simulated statistics.
        for model = 1:length(modeltypes)
            statistics(model+1,:,1) = squeeze(simulated_stats(model,1,:,t));
            statistics(model+1,:,2) = squeeze(simulated_stats(model,2,:,t));
            statistics(model+1,:,3) = squeeze(simulated_stats(model,3,:,t));
            statistics(model+1,:,5) = squeeze(simulated_stats(model,5,:,t));
            statistics(model+1,:,6) = squeeze(simulated_stats(model,6,:,t));
            statistics(model+1,:,7) = squeeze(simulated_stats(model,7,:,t));

            total_nodal_edge_length_simulated = zeros(nroi,1);
            for n = 1:nroi
                s = squeeze(simulated_stats(model,4,n,t));
                total_nodal_edge_length_simulated(n) = sum(norm(n-s));
            end
            statistics(model+1,:,4) = total_nodal_edge_length_simulated';
        end

        for model = 1:length(modeltypes)+1            
            % Now correlate the nodal statistics within each simulated model for
            % each of the 1000 simulations.
            correlated_simulated_observed_across_models(model,:,:) = corr(squeeze(statistics(model,:,:)));
        end

        % And initialise an output variable for topological dissimilarity.
        dissimilarity = zeros(length(modeltypes),1);

        for model = 1:length(modeltypes)
            % The larger the dissimilarity, the worse the model fit. Assign
            % dissimilarity to the dissimilarity_repeat field.
            dissimilarity(model,1) = norm(squeeze(correlated_simulated_observed_across_models(model+1,:,:)) - squeeze(correlated_simulated_observed_across_models(1,:,:)));
        end

        % Assign dissimilarity to the output variable.
        dissimilarity_repeat(:,t) = dissimilarity;
    end
    
    % Now test for differences between models in topological dissimilarity,
    % for this parcellation. Prepare a vector to capture model types. 
    md=zeros(1,5000);
    md(1:1000)=1;
    md(1001:2000)=2;
    md(2001:3000)=3;
    md(3001:4000)=4;
    md(4001:5000)=5;

    % Test for differences between models in topological dissimilarity.
    [p,stats,tbl]=anova1([dissimilarity_repeat(1,:), dissimilarity_repeat(2,:),...
        dissimilarity_repeat(3,:),dissimilarity_repeat(4,:),dissimilarity_repeat(5,:)],md);
    % And the multiple comparison test...
    multcompare(tbl)
    clear p stats tbl 

    %% PART 7 - COMPARE GNM RULES THROUGH SIMILARITY OF NODAL DISTRIBUTION %%
    % The final evaluative measure we use is the similarity in spatial
    % distribution of nodal degree in the empirical and simulated networks. To
    % do this, correlate nodal degree in each simulation with the observed
    % connectome. 
    
    % Correlate observed and simulated degree for each model, parcellation, and
    % run!
    for t = 1:nrep
        for model = 1:length(modeltypes)
            r = corrcoef(x{parcellation,1},squeeze(simulated_stats(model,1,:,t)));
            spatial_repeat(model,t) = r(1,2);
        end
    end

    % Now test for differences between models in capturing observed degree.
    [p,stats,tbl]=anova1([spatial_repeat(1,:) spatial_repeat(2,:) spatial_repeat(3,:) spatial_repeat(4,:) spatial_repeat(5,:)],md);
    % And the multiple comparison test...
    multcompare(tbl)
    clear p stats tbl

end

%% PART 8 - VISUALISE ENERGY LANDSCAPES - FIGURE 2D %%
% For each model, visualise the landscape. This shows that even though the
% homophily models have a larger spread of energy across simulations, the
% optimal simulations search a narrow convex optimisation parameter window.

% NOTE: Running the full 99,856 simulations for each model may take several
% hours. Therefore, you can load the following data which holds the energy
% values for all simulations. If you run your own simulations, the optimal
% parameters may vary marginally to our estimates, due to the probabilistic
% nature of the GNMs. Note that we visualise the energy landscapes for the
% Schaefer 100-node parcellation, as this was our main analysis. 
group_gnm_energy = load('data/group_gnm_energy_99856.mat');
group_structural_gnm_output.energy = group_gnm_energy.original_simulation_energy';

for model = 1:length(modeltypes)
    imagesc(reshape(group_structural_gnm_output.energy(model,:),[sqrt(nparams) sqrt(nparams)]));
    axis off
    c = colorbar('southoutside');
    c.Ticks = [0 0.5 1];
    caxis([0 1]);
    c.Label.String = 'Energy';
    c.FontSize = 20;
    c.TickLength = 0;
    % Note that we add in the parameter limits in BioRender as this results
    % in nicer formatting.
    exportgraphics(gca, sprintf('%s_energy_landscape.png',modeltypes{model}),'Resolution',600);
end

%% PART 9 - Compare Energy from Top-25 Simulations Across Models in the Schaefer-100 Node Parcellation %%
% For each model, sort the energy from lowest to highest, extract the 25
% lowest-energy simulations, and compare them across models using a one-way
% ANOVA and post-hoc tests.
sptl_sorted = sort(group_structural_gnm_output.energy(1,:));
neighbors_sorted = sort(group_structural_gnm_output.energy(2,:));
matching_sorted = sort(group_structural_gnm_output.energy(3,:));
clu_avg_sorted = sort(group_structural_gnm_output.energy(4,:));
deg_avg_sorted = sort(group_structural_gnm_output.energy(5,:));

% Create a vector that corresponds to model
md=zeros(1,125); 
md(1:25)=1;
md(26:50)=2;
md(51:75)=3;
md(76:100)=4;
md(101:125)=5;

% One-way ANOVA to compare the energies
[p,tbl,stats]=anova1([sptl_sorted(1:25) neighbors_sorted(1:25) matching_sorted(1:25) ...
    clu_avg_sorted(1:25) deg_avg_sorted(1:25)],md);
% Multiple comparisons test to assess which group differences are
% statistically significant.
multcompare(stats)

%% PART 10 - VISUALISE TOPOLOGICAL DISSIMILARITY FINGERPRINTS %%
% This corresponds to Figure 2e. Starting with the observed connectome, and
% then for each GNM model, correlate each set of measures within each
% model, after initialising an output array. Note that topological
% dissimilarity (TD) plots for each model were arranged in BioRender and
% the formatting of titles changed etc. Use the first of the 1000
% simulations for the Schaefer 100-node parcellation to provide an example 
% of the overall statistical trend. 

% Load the simulated stats for the Schaefer 100-node parcellation.
simulated_stats = load('data/simulated_statistics_lowest_energy_eta_gamma_1000_runs_schaefer-100.mat');
simulated_stats = simulated_stats.y;

% Initialise a statistics variable.
nroi = 100;
statistics = zeros(length(modeltypes)+1,total_number_of_comparison_measures,nroi);
% Assign the observed statistics.
statistics(1,1,:) = x{1,1};
statistics(1,2,:) = x{1,2};
statistics(1,3,:) = x{1,3};
statistics(1,5,:) = x{1,5};
statistics(1,6,:) = x{1,6};
statistics(1,7,:) = x{1,7};

% Get the total nodal edge length for the observed edge length.
total_nodal_edge_length_empirical = zeros(nroi,1);
for n = 1:nroi
    s = x{1,4}(n,:);
    total_nodal_edge_length_empirical(n) = sum(norm(n-s));
end
statistics(1,4,:) = total_nodal_edge_length_empirical';

% Assign the simulated statistics. Note that we have already calculated
% total nodal edge length for the simulated statistics!
for model = 1:length(modeltypes)
    for measure = 1:total_number_of_comparison_measures
        statistics(model+1,measure,:) = squeeze(simulated_stats(model,measure,:,1));
    end
end

% Now calculate the topological dissimilarity fingerprint for each model!
td_fingerprint = zeros(length(modeltypes)+1,total_number_of_comparison_measures,total_number_of_comparison_measures);
measure_list = ["Degree", "Clustering", "Betweenness", "Edge Length", "Local Efficiency", "Eigenvector", "Modularity"];
for model = 1:length(modeltypes)
    for i = 1:total_number_of_comparison_measures
        for j = 1:total_number_of_comparison_measures
            correlated_stats = corrcoef(statistics(model+1,i,:), squeeze(statistics(1,j,:)));
            td_fingerprint(model+1,i,j) = correlated_stats(1,2);
        end
    end
end

% Now correlate within the observed connectome only!
for i = 1:total_number_of_comparison_measures
    for j = 1:total_number_of_comparison_measures
        correlated_stats = corrcoef(statistics(1,i,:),statistics(1,j,:));
        td_fingerprint(1,i,j) = correlated_stats(1,2);
    end
end

% Format the different model names nicely.
formatted_modeltypes = ["Observed", "Spatial", "Neighbours", "Matching", "Clu-Avg", "Deg-Avg"];
for model = 1:length(formatted_modeltypes)
    % Visualise the fingerprint for each model.
    imagesc(squeeze(td_fingerprint(model,:,:)));
    [t,s] = title(formatted_modeltypes{model});
    t.FontSize = 20;
    t.FontWeight = 'bold';
    set(gca, 'TickLength', [0 0], 'XTickLabel', [], 'YTickLabel', ...
        measure_list, 'FontWeight', 'bold', 'FontSize', 20);
    colorbar('FontSize', 15, 'TickLength', [], 'Ticks', [-1,0,1]);
    caxis([-1 1]);
    saveas(gca,sprintf('%s_TD.png',formatted_modeltypes{model}))
end

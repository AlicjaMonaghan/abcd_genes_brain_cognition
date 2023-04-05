% This function runs the neighbour-homophily structural generative network
% model (GNM) for a single participant, using the Schaefer 100-node
% 17-network parcellation. It can be used to produce comparable results as
% results sub-section 4 in our manuscript. The function takes 5 arguments -
% participant index, the lower and upper limits of eta (distance penalty) 
% and gamma (wiring value), respectively, and the number of runs (eta-gamma
% combinations). Note that in our study, we used eta limits of -7 and -.2,
% and gamma limits of -.0667 and .6, respectively. This was a parameter
% window which bound the lowest energy group-level GNM simulations. Our
% study used 74,529 unique eta-gamma combinations for each participant.
% Correspondence to Alicja Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk

% STEPS:
% 1. Setting up the work space
% 2. Initialising GNMs
% 3. Running the grid search, collecting statistics, and calculating
% energy!

function individual_structural_neighbours_output = run_individual_neighbours_structural_gnm_grid_search(participant_index,eta_lower_limit,eta_upper_limit,gamma_lower_limit,gamma_upper_limit,nruns)
% This script details running of the neighbour homophily structural
% generative network model (GNM), using the best-performing eta and gamma 
% narrow neighbour window limits, as determined from the group consensus 
% GNM. Note that we shall be using the Schaefer-100 node 17-network
% parcellation!
%% PART 1 - Setting up the workspace.
% Change working directory. <-- SET THIS TO WHERE YOU SAVED THE DIRECTORY!
cd('abcd_genomic_variation_structural_generative_mechanisms_open/');
% Add path for the brain connectivity toolbox 
brainconnectivity_path = 'toolboxes_and_functions/2019_03_03_BCT';
addpath(brainconnectivity_path);
% Add path for the energy function i.e. computing Kolmogorov-Smirnov
% statistics for distributions of 4 graph theory metrics between observed
% and GNM-simulated connectomes. 
functions_path = 'toolboxes_and_functions/';
addpath(functions_path);
% Now load the parcellation coordinates. 
schaefer100_metadata = load('data/schaefer100x17_1mm_info.mat');
schaefer100_coordinates = [schaefer100_metadata.schaefer100x17_1mm_info.x_mni,...
    schaefer100_metadata.schaefer100x17_1mm_info.y_mni,... 
    schaefer100_metadata.schaefer100x17_1mm_info.z_mni];
% Calculate the euclidean distance between the coordinates.
D = squareform(pdist(schaefer100_coordinates));

% Load up the individual target connectomes and seed. 
target_connectomes = load('data/example_individual_gnm_targets.mat');
target_connectomes = target_connectomes.abcd_thresholded_27_streamlines;
participant_target_connectome = squeeze(target_connectomes(participant_index,:,:));
seed = load('data/example_seed_gnm.mat');
seed = seed.seed_network;
% Update user about which participant we are processing.
fprintf('Processing participant %d.\n', participant_index);

%% PART 2 - Initialise Structural Generative Network Models
disp('Setting up model parameters.\n');
modelvar = [{'powerlaw'},{'powerlaw'}]; 
% Plot 2D coordinates, spanning the eta and gamma limits, respectively. 
% This grid will represent different combinations of eta and gamma, which 
% shall be evaluated below (in terms of KS-statistics and energies). Eta
% and gamma are inputted by the user. 
[p,q] = meshgrid(linspace(eta_lower_limit,eta_upper_limit,sqrt(nruns)),...
                   linspace(gamma_lower_limit,gamma_upper_limit, ...
                   sqrt(nruns)));
params = unique([p(:) q(:)],'rows');
nparams = size(params,1);
%% PART 3 - Run the Grid Search!
% We shall run a grid search for the neighbour model for each participant,
% using the parameters defined above. In order to calculate topological
% dissimilarity, we shall keep edge-wise and node-wise statistics, such as
% K, parameterised K, and wiring probability.

% Get key observed statistics
comparison_measures_for_ks = 4;

% These LOCAL measures will be used to calculate topological dissimilarity
% Number of bi-directional connections in the target
m = nnz(participant_target_connectome)/2; 
% Network cardinality
n = length(participant_target_connectome); 
% Create a cell to hold graph theory metrics of the observed connectomes
x = cell(8,1);
% Calculate degree. 
x{1} = sum(participant_target_connectome,2); 
% Find the clustering coefficient of the connectome at a nodal level
x{2} = clustering_coef_bu(participant_target_connectome); 
% Find the betweenness centrality
x{3} = betweenness_bin(participant_target_connectome)'; 
% Find the upper triangle of the connectome larger than 0, 
% and make into a double (EDGE LENGTH)...
x{4} = D(triu(participant_target_connectome,1) > 0);
% And local efficiency...
x{5} = efficiency_bin(participant_target_connectome,1);
% And eigenvector centrality...
x{6} = eigenvector_centrality_und(participant_target_connectome);
% And nodal participant coefficient...
x{7} = participation_coef(participant_target_connectome,0);
% And modularity...
x{8} = modularity_und(participant_target_connectome);

% Initialise the output variables. 
individual_structural_neighbours_output = struct;
individual_structural_neighbours_output.energy =  zeros(nparams,1);
individual_structural_neighbours_output.KS = zeros(nparams,comparison_measures_for_ks);

% Now run the grid search with each structural GNM! Create a counter to
% monitor how many runs have been processed, and update the user at each
% point!
start_PP_loop = tic;
% Run the GNM.
B = generative_model(seed,D,m,"neighbors",modelvar,params);
disp('Starting to run the generative neighbours model');
% Find the number of runs.
nB = size(B,2);
% Create a matrix of zeros for each measure of model fit
K = zeros(nB,comparison_measures_for_ks); 
% For the number of GNM parameter combinations... 
for iB = 1:nB 
 b = zeros(n); 
 b(B(:,iB)) = 1; 
 b = b + b'; 

 % Create a cell to hold the synthetic simulated statistics
 % needed to calculate the energy!
 y = cell(comparison_measures_for_ks,1);

 % Note that we will only collect the statistics needed to calculate the KS
 % statistics because for the best-fitting model at the best-fitting 
 % parameter combination, we will re-run the model with that specific 
 % combination, and calculate these statistics. We have calculated a full 
 % set of 8 nodal statistics for the observed connectomes, because this
 % will not depend on the model.
 
 % Sum the number of connections in each column (number of streamlines), 
 % hence DEGREE
 y{1} = sum(b,2);
 % Find clustering coefficient...
 y{2} = clustering_coef_bu(b); 
 % And betweenness-centrality...
 y{3} = betweenness_bin(b)'; 
 % Find the upper triangle of the connectome larger than 0, and make into 
 % a double (EDGE LENGTH)
 y{4} = D(triu(b,1) > 0); 
 
 % Calculate the energy of the synthetic network, using degree, clustering
 % coefficient, betweenness centrality, and edge length distributions 
 % (see Betzel et al., 2016).
 for j = 1:comparison_measures_for_ks 
    K(iB,j) = fcn_ks(x{j},y{j}); 
 end
individual_structural_neighbours_output.KS = K; % Keep ks statistics
% This is the energy of the model - the lower the energy, the better the 
% model fit. Note that the distribution of energy values does not matter,
% only the smallest energy value.
individual_structural_neighbours_output.Energy = max(K,[],2); 
end

% Now update the user about how long each model took overall!    
model_processing_complete = toc(start_PP_loop); 
% Time how long each model has taken to process
fprintf('Completed the neighbours model for participant %d in %g seconds!\n', ...
    participant_index,round(model_processing_complete,2));
% And save in the data directory
save(sprintf('data/individual_gnm_neighbours_participant_%d.mat',participant_index),"individual_structural_neighbours_output");
end

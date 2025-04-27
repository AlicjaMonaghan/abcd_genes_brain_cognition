% To probe the effect of changing parameters on stochasticity and network
% architecture, we conduct a simulation experiment where we simulate
% connectivity at decreasing deciles of the optimal group-level eta and
% gamma estimates. These analyses are inspired by Sofia Carozza's 2023 
% Developmental Psychobiology paper "Early adversity changes the economic
% conditions of mouse structural brain organization". Written by Alicja
% Monaghan in February 2025.

%% PART 1 - Set up the workspace
clear;clc;
cd('/Users/alicjamonaghan/Desktop/abcd_genes_brain_cognition/');
% Add path to the brain connectivity toolbox
addpath('2019_03_03_BCT/');
% Load the group-level seed. The same seed is used for the group-level and
% PGS contrast analyses.
seed = load('data/seed_across_parcellations.mat').seed.schaefer100;
% Find the number of bidirectional connections in the seed
mseed = nnz(seed)/2;
% Load the coordinates for the Schaefer 100-node parcellation, and
% calculate the Euclidean distance between points. 
schaefer100x17_1mm_info = load('data/schaefer100x17_1mm_info.mat');
schaefer100x17_1mm_info = schaefer100x17_1mm_info.schaefer100x17_1mm_info;
parcellation_coordinates = [schaefer100x17_1mm_info.x_mni,...
            schaefer100x17_1mm_info.y_mni, schaefer100x17_1mm_info.z_mni];
D = squareform(pdist(parcellation_coordinates));
% Load the target we're aiming to simulate, and get the number of
% bi-directional connections.
target = load('data/group_target_across_parcellations.mat').target.schaefer100;
m = nnz(target)/2;
% Set the minimum edge value
epsilon = 1e-5;
% And the number of regions
n = length(D);
%% Part 2 - Examine properties of increasingly stochastic networks
% Across deciles of the best-fitting group-level eta and gamma parameters,
% calculate the global efficiency and quantify randomness by topological
% dissimilarity to randomly rewired connectomes.
eta_params = [-3.0610, -2.9630];
gamma_params = [.2138, .2144];
deciles = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
nruns = 1000;
synthetic_connectome_global_efficiency = zeros(length(deciles), nruns);
% We're using the optimal eta and gamma parameters at the group level, with
% the same seed and target as before.
eta_deciles = deciles*-2.911;
gamma_deciles = deciles*.244;
% Compute topological dissimilarity of simulations to randomly
% degree-preserved networks
rand_dissim_pct = zeros(length(deciles), nruns);
for decile_idx = 1:length(deciles)
    for run = 1:nruns
        B = generative_model(seed, D, m, 'neighbors', {'powerlaw', 'powerlaw'}, ...
            [eta_deciles(decile_idx) gamma_deciles(decile_idx)]);
        b = zeros(n);
        b(B(:, 1)) = 1;
        b = b + b';
        % Calculate global efficiency of synthetic connectomes
        synthetic_connectome_global_efficiency(decile_idx, run) = efficiency_bin(b);
        % For topological dissimilarity, calculate degree, clustering,
        % betweenness-centrality, edge length, local efficiency and
        % eigenvector centrality.
        nodedeg = degrees_und(b)';
        nodeclust = clustering_coef_bu(b);
        nodebc= betweenness_bin(b)';
        r = triu(b,1)>0;
        total_edge_length = zeros(n, 1);
        for roi=1:n
            s = r(roi, :);
            total_edge_length(roi) = sum(D(roi,s));
        end
        nodelocalefficiency = efficiency_bin(b,1);
        nodeec = eigenvector_centrality_und(b);
        tostats_synth = cat(2,nodedeg,nodeclust,nodebc,total_edge_length,nodelocalefficiency,nodeec);
        tomatrix_synth= corr(tostats_synth);
        clear b B nodedeg nodeclust nodebc r total_edge_length s nodelocalefficiency nodeec
        % Re-wire each edge once from the empirical connectome, whilst 
        % preserving the degree distribution
        connectome_rand = randmio_und(target, 1);
        % Calculate the same attributes as above
        nodedeg = degrees_und(connectome_rand)';
        nodeclust = clustering_coef_bu(connectome_rand);
        nodebc = betweenness_bin(connectome_rand)';
        r = triu(connectome_rand,1)>0;
        total_edge_length = zeros(n, 1);
        for roi=1:n
            s = r(roi, :);
            total_edge_length(roi) = sum(D(roi,s));
        end
        nodelocalefficiency = efficiency_bin(connectome_rand,1);
        nodeec = eigenvector_centrality_und(connectome_rand);
        % Concatenate the arrays and correlate
        tostats_rand = cat(2,nodedeg,nodeclust,nodebc,total_edge_length,nodelocalefficiency,nodeec);
        tomatrix_rand = corr(tostats_rand);
        % And calculate the topological dissimilarity
        rand_dissim_pct(decile_idx, run) = norm(tomatrix_synth-tomatrix_rand);
        fprintf('Run %g of %g for decile %g\n', run, nruns, decile_idx);
    end
end
% Save the synthetic connectome global efficiency and the topological
% dissimilarity compared to randomly-rewired networks
save('data/synthetic_connectome_decile_global_efficiency.mat', "synthetic_connectome_global_efficiency")
save('data/topological_dissimilarity_random_rewiring.mat', 'rand_dissim_pct');

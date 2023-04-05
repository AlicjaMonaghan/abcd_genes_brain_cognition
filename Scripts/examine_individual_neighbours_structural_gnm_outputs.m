% This script details pulling output for the neighbours structural 
% generative network model for each participant, using eta and gamma limits
% determined by the best-fitting model from a group model run across all 
% 2153 ABCD participants with a consensus-based distance-dependent network. 
% We shall sort the model energy from lowest to highest, and extract the
% best-fitting model parameters for each participant. We shall re-run the 
% neighbour structural GNM for the best-fitting parameter combination for 
% each participant, in order to collect statistics throughout each run, 
% allowing us to plot parameterised nodal costs and values. The
% individual-level outputs for the neighbours structural GNM were obtained
% by running run_individual_neighbours_structural_gnm_grid_search.mat, with
% eta limits of -7 and -.20, and gamma limits of -.0667 and .60,
% respectively, with 75,000 runs (corresponding to 74,529 unique eta-gamma
% combinations). Correspondence should be directed to Alicja Monaghan,
% alicja.monaghan@mrc-cbu.cam.ac.uk

% STEPS:
% 1. Setting up the workspace.
% 2. Visualising optimal individual-level eta-gamma combinations on the
% energy landscape - FIGURE 3A. 
% 3. Re-running lowest-energy individual simulations to extract simulated
% nodal statistics to plot using the
% abcd_group_and_individual_gnm_visualisation.R script.
%% PART 1 - Set up the workspace.
% Set the current working directory. <-- SET THIS TO WHERE YOU SAVED THE
% DIRECTORY!
clear; clc;
cd('/abcd_genomic_variation_structural_generative_mechanisms_open/'); 
% Add path for the brain connectivity toolbox 
brainconnectivity_path = 'toolboxes_and_functions/2019_03_03_BCT';
addpath(brainconnectivity_path);
% Plot 2D coordinates, spanning the eta and gamma limits, respectively. 
% This grid will represent different combinations of eta and gamma, which 
% shall be evaluated below (in terms of KS-statistics and energy).
etalimits = [-7,-.2000];
gamlimits = [-0.0667,.6000];
nruns = 75000;
[p,q] = meshgrid(linspace(etalimits(1),etalimits(2),sqrt(nruns)),...
                   linspace(gamlimits(1),gamlimits(2),sqrt(nruns)));
params = unique([p(:) q(:)],'rows');
nparams = size(params,1);
% Load the energy values for all 74,529 unique eta-gamma combinations for
% each of 2153 participants. Note that whilst we initially used 2154 
% participants, we removed one participant whose mean connectome density
% (.02%) was considerably smaller than that of the remaining participants 
% (6.55%), suggesting incorrect connectome reconstruction.
individual_gnm_energy = load('data/individual_gnm_energy.mat');
individual_gnm_energy = individual_gnm_energy.individual_gnm_energy;
nsub = size(individual_gnm_energy,1);
% Load the Schaefer 100-node 17-network parcellation meta-data.
schaefer100_metadata = load('data/schaefer100x17_1mm_info.mat');
% Now load the parcellation coordinates. 
schaefer100_coordinates = [schaefer100_metadata.schaefer100x17_1mm_info.x_mni,...
    schaefer100_metadata.schaefer100x17_1mm_info.y_mni,... 
    schaefer100_metadata.schaefer100x17_1mm_info.z_mni];
% Calculate the euclidean distance between the coordinates.
D = squareform(pdist(schaefer100_coordinates));
% And find the number of regions in the parcellation.
nroi = 100;
% Load the GNM target connectomes for all participants.
target_connectomes = load('data/individual_gnm_targets.mat');
target_connectomes = target_connectomes.abcd_thresholded_27_streamlines;
% Load the seed.
seed = load('data/example_seed_gnm.mat');
seed = seed.seed_network;

%% PART 2 - Visualise Optimal Eta-Gamma Combinations on Energy Landscape %%
% This corresponds to Figure 3a, which is a scatter plot of individual
% eta-gamma combinations producing the lowest-energy (best-fitting) GNM.
% Fist, for each participant, find the index of their lowest energy 
% simulation, and extract the corresponding eta and gamma parameters.
optimal_eta = zeros(nsub,1);
optimal_gamma = zeros(nsub,1);
lowest_energy = zeros(nsub,1);
for sub = 1:nsub
    [M,I] = min(individual_gnm_energy(sub,:));
    lowest_energy(sub,1) = M;
    optimal_eta(sub,1) = params(I,1);
    optimal_gamma(sub,1) = params(I,2);
end
% Now plot the energy landscape i.e. individual energy values and the eta
% and gamma limits they correspond to.
imagesc(squeeze(mean(reshape(individual_gnm_energy, [nsub sqrt(nparams) sqrt(nparams)]),1)));
xlabel('\fontsize{30} \eta','Interpreter','tex', 'FontWeight','bold', 'Position',[150 270])
ylabel('\fontsize{30} \gamma','Interpreter','tex', 'FontWeight','bold','Position',[-30 150], 'Rotation', 360)
xticks([2 270]); 
xticklabels({round(min(optimal_eta(:)),3), round(max(optimal_gamma(:)),3)});
set(gca, 'TickLength', [0 0])
a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'FontSize', 15)
yticks([2 270]); 
yticklabels({round(min(optimal_gamma(:)),3), round(max(optimal_gamma(:)),3)});
axis square
box off
hold on
c = colorbar;
set(c, 'TickLength', 0, 'FontSize', 15)
set(get(c, 'ylabel'), 'string', 'Energy', 'FontSize', 20, 'FontWeight', 'bold', 'Rotation', 90)
hold off
% Rescale the best individual eta and gamma values to the limits of the
% energy landscape plot. The x and y limits are the same. To create Figure
% 3a, we superimposed the scatter plot below onto the energy landscape
% above.
x1 = xlim;
individual_optimal_parameters_scaled = horzcat(rescale(optimal_eta, x1(1), x1(2)), ...
    rescale(optimal_gamma, x1(1), x1(2)));
scatter(individual_optimal_parameters_scaled(:,1), individual_optimal_parameters_scaled(:,2), ...
    80, 'x', 'MarkerEdgeColor', 'black', 'LineWidth', 1.5);
axis off

%% PART 3 - Re-Run Lowest-Energy GNM Simulations for Each Participant %%
% Re-run the neighbours-homophily GNM for each participant with their
% lowest-energy eta and gamma combinations. We collect the parameterised
% nodal wiring costs and values for each iteration in the simulation. We
% also collect the following graph theory metrics: degree, edge length, 
% betweenness-centrality, clustering-coefficient, local efficiency, 
% modularity, eigenvector centrality, and participation coefficients. 
% First, initialise output arrays. 
lowest_energy = zeros(nsub,1);
nodal_parameterised_wiring_value = zeros(nsub,nroi,nroi);
nodal_parameterised_wiring_cost = zeros(nsub,nroi,nroi);
% Note that the first 'slice' of the fourth dimension of this array
% indicates the observed connectome, and the second 'slice' indicates the
% simulated connectome.
lowest_energy_grapth_theory_metrics = zeros(nsub,nroi,8,2);
% Set the number of comparison measures needed to compute energy.
comparison_measures_for_ks = 4;
% Set some model specifications.
modelvar = [{'powerlaw'},{'powerlaw'}]; 
epsilon = 1e-5;
% Re-run the GNM!
for sub = 1:20
    fprintf('Collecting parameterised statistics for participant %d.\n',sub);
    % Initialise output variables. K is nodal wiring value, Fk is 
    % parameterised nodal wiring. P is wiring probability, and A are the
    % networks at each iteration. 
    Kall        = [];
    Fkall       = [];
    Pall        = [];
    Aall        = [];
    % Extract the target connectome for this participant.
    participant_target_connectome = squeeze(target_connectomes(sub,:,:));
    % Number of bi-directional connections in the target
    m = nnz(participant_target_connectome)/2; 
    % Network cardinality
    n = length(participant_target_connectome); 
    % Set the seed as A - this is the first network.
    A = seed;
    % Set K, as per Brain Connectivity Toolbox notation...
    K = (A*A).*~eye(n);
    % Keep the first K and network A
    Kall(1,:,:) = K;
    Aall(1,:,:) = A;
    % Now update K with epsilon, denoting the minimum likelihood of wiring
    K = K + epsilon;
    % Set additional parameters
    n = length(D);
    mseed = nnz(A)/2;  
    % Include all non-zero connections in the seed
    A = A > 0;
    % Switch model variables depending on whether we selected an
    % exponential or power law.
    mv1 = modelvar{1};
    mv2 = modelvar{2};
    switch mv1
        case 'powerlaw'
            Fd = D.^optimal_eta(sub,1);
        case 'exponential'
            Fd = exp(optimal_eta(sub,1)*D);
    end
    switch mv2
        case 'powerlaw'
            %gam = abs(gam);
            Fk = K.^optimal_gamma(sub,1);
        case 'exponential'
            Fk = exp(optimal_gamma(sub,1)*K);
    end
    Ff = Fd.*Fk.*~A;
    [u,v] = find(triu(ones(n),1));
    indx = (v - 1)*n + u;
    P = Ff(indx);
    % Save the first parameterised K in Fk
    Fkall(1,:,:)  = Fk;
    % Save the first probabilities
    Ff(isinf(Ff)) = 0;
    Pall(1,:,:)   = Ff;
    % Update the step...
    step = 2; 
    for i = (mseed + 1):m
        C = [0; cumsum(P)];
        r = sum(rand*C(end) >= C);
        uu = u(r);
        vv = v(r);
        x = A(uu,:);
        y = A(:,vv);
        A(uu,vv) = 1;
        A(vv,uu) = 1;
        K(uu,y) = K(uu,y) + 1;
        K(y,uu) = K(y,uu) + 1;
        K(vv,x) = K(vv,x) + 1;
        K(x,vv) = K(x,vv) + 1;
        switch mv2
            case 'powerlaw'
                Ff(uu,y) = Fd(uu,y).*(K(uu,y).^optimal_gamma(sub,1));
                Ff(y,uu) = Ff(uu,y)';
                Ff(vv,x) = Fd(vv,x).*(K(vv,x).^optimal_gamma(sub,1));
                Ff(x,vv) = Ff(vv,x)';
            case 'exponential'
                Ff(uu,y) = Fd(uu,y).*exp(K(uu,y)*optimal_gamma(sub,1));
                Ff(y,uu) = Ff(uu,y)';
                Ff(vv,x) = Fd(vv,x).*exp(K(vv,x)*optimal_gamma(sub,1));
                Ff(x,vv) = Ff(vv,x)';
        end
        Ff(A) = 0;
        P = Ff(indx);
        Kall(step,:,:)  = K;                   
        Fkall(step,:,:) = Fk;                  
        Aall(step,:,:)  = A;    
        Pall(step,:,:)  = Ff;                
        % Change the step
        step = step+1;
    end
    % Assign the parameterised nodal wiring value to the output variable.
    nodal_parameterised_wiring_value(sub,:,:) = squeeze(mean(Fkall,1));
    % And the parameterised nodal wiring cost...
    Fd(isinf(Fd)) = 0;
    nodal_parameterised_wiring_cost(sub,:,:) = Fd;
    % B is the actual generative model output!
    B = find(triu(A,1));
    % Create a cell array to hold graph theory metrics for the observed and 
    % simulated networks.
    x = cell(comparison_measures_for_ks,1);
    % Calculate degree and assign to relevant output arrays. 
    x{1} = sum(participant_target_connectome,2); 
    lowest_energy_grapth_theory_metrics(sub,:,1,1) = x{1};
    % Find the clustering coefficient of the connectome at a nodal level
    x{2} = clustering_coef_bu(participant_target_connectome); 
    lowest_energy_grapth_theory_metrics(sub,:,2,1) = x{2};
    % Find the betweenness centrality
    x{3} = betweenness_bin(participant_target_connectome)'; 
    lowest_energy_grapth_theory_metrics(sub,:,3,1) = x{3};
    % Find nodal edge length.
    nodal_total_edge_length = sum(D.*(participant_target_connectome>0));
    for n = 1:nroi
        s = (n);
        x{4}(1,n) = sum(norm(n-s));
    end
    % And local efficiency...
    lowest_energy_grapth_theory_metrics(sub,:,5,1) = efficiency_bin(participant_target_connectome,1);
    % And eigenvector centrality...
    lowest_energy_grapth_theory_metrics(sub,:,6,1) = eigenvector_centrality_und(participant_target_connectome);
    % And nodal participant coefficient...
    lowest_energy_grapth_theory_metrics(sub,:,7,1) = participation_coef(participant_target_connectome,0);
    % And optimal community structure to maximise modularity...
    lowest_energy_grapth_theory_metrics(sub,:,8,1) = modularity_und(participant_target_connectome);
    % Create a matrix of zeros for the cardinality of the network.
    b = zeros(n); 
    % For iB indexing the parameter combination in the GNM, find the 
    % associated column in the cardinality network. 
    b(B(:,1)) = 1;
    % Represent b as an adjacency matrix 
    b = b + b'; 
    % Initialise cell for measures needed to calculate energy.
    y = cell(comparison_measures_for_ks,1);
    % Sum the number of connections in each column (number of
    % streamlines). This is degree.
    y{1} = sum(b,2); 
    lowest_energy_grapth_theory_metrics(sub,:,1,2) = y{1};
    % Clustering coefficient.
    y{2} = clustering_coef_bu(b); 
    lowest_energy_grapth_theory_metrics(sub,:,2,2) = y{1};
    % Betweenness-centrality.
    y{3} = betweenness_bin(b)'; 
    lowest_energy_grapth_theory_metrics(sub,:,3,2) = y{1};
    % Nodal edge length.
    y{4} = D(triu(b,1) > 0);
    % To plot nodal edge length, we need to find the TOTAL edge length from
    % each node.
    total_node_edge_length = zeros(nroi,1);
    for n = 1:nroi
        s = y{4}(n,1);
        total_node_edge_length(n,1) = sum(norm(n-s));
    end
    y{4} = total_node_edge_length';
    % Local efficiency...
    lowest_energy_grapth_theory_metrics(sub,:,5,2) = efficiency_bin(b,1);
    % Eigenvector centrality...
    lowest_energy_grapth_theory_metrics(sub,:,6,2) = eigenvector_centrality_und(b);
    % Nodal participant coefficient...
    lowest_energy_grapth_theory_metrics(sub,:,7,2) = participation_coef(b,0);
    % Optimal community structure to maximise modularity...
    lowest_energy_grapth_theory_metrics(sub,:,8,2) = modularity_und(b);
end

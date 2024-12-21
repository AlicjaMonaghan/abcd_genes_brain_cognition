% This script details pulling output from structural GNMs for the consensus
% group network, representing/averaged across all 2154 ABCD participants,
% determining the best-fitting parameters across for the group. We shall 
% first collect energies, and sort them per model per parcellation. We 
% shall then use the associated parameters as the boundary for the narrow
% matching window, for which the best-fitting model shall be run. To
% calculate topological similarity, we shall pull the synthetic output for
% each model and parcellation, and calculate local and global statistics
% per run. 

%% PART 1 - Set up the workspace.
% Load up the directory containing the structural GNM outputs.
clear; clc;
Structural_GNM_Outputs_RootDir = '\\cbsu\data\Imaging\projects\external\abcd\analyses\Alicja\Structural_GNMs_Alicja_V3\Group_GNMs'
addpath(Structural_GNM_Outputs_RootDir)
cd(Structural_GNM_Outputs_RootDir);
% Add the path for the Brain Connectivity Toolbox
brainconnectivity_path = '/imaging/astle/am10/Toolboxes/2019_03_03_BCT';
addpath(brainconnectivity_path);
% Specify the number of runs we'd like to extract data for.
nruns = 100000;
% Load up the different model names and specify the number tested 
modeltypes = ["sptl", "neighbors", "matching", "clu-avg", "deg-avg"];
nmodels = length(modeltypes);
% Load up the grid search parameter limits
etalimits = [-7,7];
gamlimits = [-7,7];
[p,q]     = meshgrid(linspace(etalimits(1),etalimits(2),sqrt(nruns)),...
                   linspace(gamlimits(1),gamlimits(2),sqrt(nruns)));
params = unique([p(:) q(:)],'rows');
nparams = size(params,1);

% Load in the different parcellation options. 
Parcellation_Options = ["schaefer100x17","brainnetome246","schaefer400x17"];
% And the measures we will be using for calculating topological
% dissimilarity (TD)
Measures = {'Degree', 'Clustering', 'Betweenness Centrality', 'Edge Length', 'Local Efficiency', 'Centrality', 'Participation Coefficient', 'Modularity'};

%% PART 2 - Pull and Sort Model Energies %%
% Due to the size of the MATLAB files, we pulled the lowest energy for each
% model for each parcellation, and their corresponding eta and gamma
% parameter, alongside the lowest-energy 1000 simulations. Therefore, we
% shall load these outputs now. Each has the same two first dimensions:
% parcellation (out of 3) x model (out of 5), using the same order as
% above.
cd('Group_Sorted_Outputs\')
load('ABCD_Minimum_Energy_Group.mat')
load("ABCD_Top_1000_Simulations_Group.mat")
load("ABCD_Best_Eta_Group.mat")
load("ABCD_Best_Gamma_Group.mat")
% We find that the lowest energy achieved for the Schaefer 400-node
% parcellation was achieved by the neighbours model, with eta and gamma
% parameters of -2.911 and .3778, respectively. Now, let's find the
% parameters which are associated with the neighbours model, for the top
% 10% lowest-energy simulations. These parameters shall be used as the
% boundaries for running the individual generative network models.
load("ABCD_Top_1000_Eta_Group.mat")
load("ABCD_Top_1000_Gamma_Group.mat")
% Eta lower boundary of -7
Eta_Lower_Boundary = min(ABCD_Top_1000_Eta(3,2,:))
% Eta upper boundary of -1.9778
Eta_Upper_Boundary = max(ABCD_Top_1000_Eta(3,2,:))
% Gamma lower boundary of -0.2889
Gamma_Lower_Boundary = min(ABCD_Top_1000_Gamma(3,2,:))
% Gamma upper boundary of 2.5556
Gamma_Upper_Boundary = max(ABCD_Top_1000_Gamma(3,2,:))

%% PART 3 - Calculating Topological Dissimilarity! %%
% As an additional check for the best-fitting model, we will examine the
% topological similarity for the best-fitting parameter combination for
% each model. We shall find the indices of the best eta and gamma
% combinations for each model and parcellation, and then select the
% corresponding synthetic network, from which we shall calculate simulated
% statistics. Since the outputs are so large for the finer parcellations,
% we shall re-run each generative model for the best-fitting parameter
% combinations.

% First, intialise a structure of cells for the outputs. The first will 
% hold the empirical statistics (number of parcellations x number of 
% measures), and the second simulated for all models and parcellations 
% (number of parcellations x number of models x number of measures).
Empirical_Statistics = cell(length(Parcellation_Options),10);
Simulated_Statistics= cell(length(Parcellation_Options),nmodels,10)
cd('/imaging/projects/external/abcd/analyses/Alicja/Structural_GNMs_Alicja_V3/Group_GNMs/')
Correlated_Simulated_Observed_Across_Models_and_Parcellations = zeros(length(Parcellation_Options),(1+length(modeltypes)),length(Measures),length(Measures));
Dissimilarity_Across_Models_and_Parcellations = zeros(length(Parcellation_Options),length(modeltypes),1);
Spatial_Embedding_Across_Models_and_Parcellations = zeros(length(Parcellation_Options),nmodels,1);

% Now run the loop, where we shall re-run all GNMs for each model and
% parcellation based on the optimal parameter combination, to retrieve
% simulated statistics. 
for Parcellation = 1
    % Load up the correct parcellation coordinates.
        if Parcellation == 1
            Coordinates = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x17_1mm_info.mat');
            D = squareform(pdist([Coordinates.schaefer100x17_1mm_info.x_mni Coordinates.schaefer100x17_1mm_info.y_mni Coordinates.schaefer100x17_1mm_info.z_mni]));
            nroi = length(D);
        elseif Parcellation == 2
            Coordinates = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/brainnetome246_info.mat');
            D = squareform(pdist([Coordinates.brainnetome246.x_mni Coordinates.brainnetome246.y_mni Coordinates.brainnetome246.z_mni]));
            nroi = length(D);
        else 
            Coordinates = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer400x17_1mm_info.mat');
            D = squareform(pdist([Coordinates.schaefer400x17_1mm_info.x_mni Coordinates.schaefer400x17_1mm_info.y_mni Coordinates.schaefer400x17_1mm_info.z_mni]));
            nroi = length(D);
        end

    % Now, to calculate topological dissimilarity, we need to load up the 
    % observed connectome (before any modelling) for each parcellation and 
    % model. Specifically, we need to calculate additional statistics, so
    % that the number of empirical and synthetic statistics match.
    Target = load(sprintf('ABCD_Group_Consensus_Based_Structural_Connectome_%s_Parcellation.mat',Parcellation_Options(Parcellation)));
    % Update user!
    fprintf('Calculating statistics for the empirical connectome for the %s parcellation.\n',Parcellation_Options(Parcellation))
    % Calculate the same statistics as above, and append to the output
    % variable!
    Empirical_Statistics{Parcellation,1} = sum(Target.Group_Binarised_Structural_Connectome,2); 
    Empirical_Statistics{Parcellation,2} = clustering_coef_bu(Target.Group_Binarised_Structural_Connectome); 
    Empirical_Statistics{Parcellation,3} = betweenness_bin(Target.Group_Binarised_Structural_Connectome)'; 
    Empirical_Statistics{Parcellation,4} = D(triu(Target.Group_Binarised_Structural_Connectome,1) > 0);
    Empirical_Statistics{Parcellation,5} = efficiency_bin(Target.Group_Binarised_Structural_Connectome,1);
    Empirical_Statistics{Parcellation,6} = eigenvector_centrality_und(Target.Group_Binarised_Structural_Connectome);
    Empirical_Statistics{Parcellation,7} = participation_coef(Target.Group_Binarised_Structural_Connectome,0);
    Empirical_Statistics{Parcellation,8} = rich_club_bu(Target.Group_Binarised_Structural_Connectome); 
    Empirical_Statistics{Parcellation,9} = efficiency_bin(Target.Group_Binarised_Structural_Connectome);
    Empirical_Statistics{Parcellation,10} = modularity_und(Target.Group_Binarised_Structural_Connectome);

    for Model = 1:nmodels

        % Update user about which parcellation and model combination are
        % being processed.
        fprintf('Running GNMs for the %s parcellation!\n',Parcellation_Options(Parcellation));
        % Load the correct seed and target networks!
        Seed = load(sprintf('ABCD_Seed_Network_%s_Parcellation.mat',Parcellation_Options(Parcellation)));
        %Target = load(sprintf('ABCD_Group_Consensus_Based_Structural_Connectome_%s_Parcellation.mat',Parcellation_Options(Parcellation)));
        % Set up the generative network models
        modelvar = [{'powerlaw'},{'powerlaw'}]; 
        params = [Best_Eta(Parcellation,Model) Best_Gamma(Parcellation,Model)];
        % Number of bi-directional connections in the target
        m = nnz(Target.Group_Binarised_Structural_Connectome)/2; 
        % Network cardinality
        n = length(Target.Group_Binarised_Structural_Connectome); 
        % Here, we are running the actual generative model!
        B = generative_model(Seed.ABCD_Seed_Network,D,m,modeltypes(Model),modelvar,params);
        fprintf('Starting to run the generative %s model for the %s parcellation!\n', modeltypes(Model), Parcellation_Options(Parcellation));
        
        % Now calculate the synthetic network measures. These LOCAL 
        % measures will be used to calculate topological dissimilarity.
        % Append these outputs to the cell structure we initialised
        % earlier.

        nB = size(B,2); %Finding the size of the second dimension of the generative model
        for iB = 1:nB %For the number of columns/parameters in the generative network model
            b = zeros(n); %Create a matrix of zeros for the cardinality of the network (size of Connectivity Dimensions)
            b(B(:,iB)) = 1; %For iB indexing the parameter of the generative network model, find the associated column in the cardinality network (size 246 x 246)
            b = b + b'; %Represent b as an adjacency matrix 
    
            % Sum the number of connections in each column (number of streamlines), hence DEGREE
            Simulated_Statistics{Parcellation,Model,1} = sum(b,2); 
            % Find the clustering coefficient of the connectome at a nodal level
            Simulated_Statistics{Parcellation,Model,2} = clustering_coef_bu(b); 
            % Find the betweenness centrality
            Simulated_Statistics{Parcellation,Model,3} = betweenness_bin(b)'; 
            % Find the upper triangle of the connectome larger than 0, and make into a double (EDGE LENGTH)
            Simulated_Statistics{Parcellation,Model,4} = D(triu(b,1) > 0);
            % Find the local efficiency
            Simulated_Statistics{Parcellation,Model,5} = efficiency_bin(b,1);
            % Find the eigenvector centrality
            Simulated_Statistics{Parcellation,Model,6} = eigenvector_centrality_und(b);
            % Finally, find the nodal participant coefficient
            Simulated_Statistics{Parcellation,Model,7} = participation_coef(b,0);
    
            % The below 3 measures are additional measures, at a GLOBAL level.
            % Rich-club coefficient (fraction of edges connecting nodes with degree k from the maximum number of possible edges a node shares)
            Simulated_Statistics{Parcellation,Model,8} = rich_club_bu(b); 
            % Calculates the global efficiency of the network
            Simulated_Statistics{Parcellation,Model,9} = efficiency_bin(b);
            % Find the modularity of the network using Newman's spectral community detection algorithm
            Simulated_Statistics{Parcellation,Model,10} = modularity_und(b);
       end
       
    end

    % Now we need to calculate topological dissimilarity (TD)! Initialise 
    % an output variable for each parcellation, where the first entry of 
    % the first dimension will be the empirical statistics, and the 
    % remainder will be the simulated!
    Statistics = zeros((1+nmodels),nroi,length(Measures));
    % Now add the empirical connectome properties for each parcellation.
    % Note that the order of the measures differs in the statistics and
    % empirical statistics data frames, because we are not including all
    % empirical statistics in the topological dissimilarity, but only those
    % with nodal values.
    Statistics(1,:,1) = Empirical_Statistics{Parcellation,1};
    Statistics(1,:,2) = Empirical_Statistics{Parcellation,2};
    Statistics(1,:,3) = Empirical_Statistics{Parcellation,3};
    Statistics(1,:,5) = Empirical_Statistics{Parcellation,5};
    Statistics(1,:,6) = Empirical_Statistics{Parcellation,6};
    Statistics(1,:,7) = Empirical_Statistics{Parcellation,7};
    Statistics(1,:,8) = Empirical_Statistics{Parcellation,10};
    % Calculate the total distance from each node
    Node_Edge_Length_Empirical = size(nroi);
    for n = 1:nroi
        s = Empirical_Statistics{Parcellation,4}(n,:);
        Node_Edge_Length_Empirical(n) = sum(norm(n-s));
    end
    Statistics(1,:,4) = Node_Edge_Length_Empirical.';

    % And now load all of the simulated model characteristics
    for model = 1:length(modeltypes)
        Statistics(model+1,:,1) = Simulated_Statistics{Parcellation,model,1};
        Statistics(model+1,:,2) = Simulated_Statistics{Parcellation,model,2};
        Statistics(model+1,:,3) = Simulated_Statistics{Parcellation,model,3};
        Statistics(model+1,:,5) = Simulated_Statistics{Parcellation,model,5};
        Statistics(model+1,:,6) = Simulated_Statistics{Parcellation,model,6};
        Statistics(model+1,:,7) = Simulated_Statistics{Parcellation,model,7};
        Statistics(model+1,:,8) = Simulated_Statistics{Parcellation,model,10};        
        % Now calculate the total distance from each node
        Node_Edge_Length_Simulated = nroi;
        for n = 1:nroi
            s = Simulated_Statistics{Parcellation,model,4}(n,:);
            Node_Edge_Length_Simulated(n) = sum(norm(n-s));
        end
        Statistics(model+1,:,4) = Node_Edge_Length_Simulated.';
        clear Node_Edge_Length_Simulated
    end

    % Now correlate the nodal statistics within each model, after 
    % initialising output arrays. The final two dimensions provides the 
    % topological dissimilarity for each model.
    Correlated_Simulated_Observed_Across_Models = zeros((1+length(modeltypes)),length(Measures),length(Measures));
    for model = 1:(1+length(modeltypes))
        Correlated_Simulated_Observed_Across_Models(model,:,:) = corr(squeeze(Statistics(model,:,:)));
    end

    % Assign to the output variable!
    Correlated_Simulated_Observed_Across_Models_and_Parcellations(Parcellation,:,:,:) = Correlated_Simulated_Observed_Across_Models;
    
    % Finally, calculate the topological dissimilarity, where we correlate the 
    % dissimilarity between each model's topological fingerprint and the 
    % empirical connectome topological fingerprint. The model with the
    % smallest topological dissimilarity is deemed as accounting for the local
    % statistics in the best way!
    Dissimilarity = zeros(length(modeltypes),1);
    for model = 1:length(modeltypes)
        Dissimilarity(model,1) = norm(squeeze(Correlated_Simulated_Observed_Across_Models(model+1,:,:)) - squeeze(Correlated_Simulated_Observed_Across_Models(1,:,:)));
    end
    fprintf('Calculated topological dissimilarity for %s parcellation!\n',Parcellation_Options(Parcellation));

    % And assign to the output variable!
    Dissimilarity_Across_Models_and_Parcellations(Parcellation,:) = Dissimilarity;

    %% PART 4 - Spatial Embedding %%
    % As an additional criterion, we assessed spatial embedding of networks,
    % namely correlations between nodal degree in the empirical and simulated
    % connectomes.
    Spatial_Embedding = zeros(nmodels,1);
    for model = 1:length(modeltypes)
        Correlation_4_by_4_Matrix = corrcoef(Statistics(1,:,1),Statistics((model+1),:,1));
        % Select off-diagonal coefficient i.e. single coefficient for a whole
        % matrix!
        Spatial_Embedding(model,1) = Correlation_4_by_4_Matrix(1,2);
    end
    fprintf('Calculated spatial embedding for %s parcellation!\n',Parcellation_Options(Parcellation));

    % And assign to the output variable
    Spatial_Embedding_Across_Models_and_Parcellations(Parcellation,:) = Spatial_Embedding;
   
end
     
%% PART 5 - Save Outputs!
save('Group_Sorted_Outputs/ABCD_Group_Simulated_Statistics_Across_Models_and_Parcellations.mat','Simulated_Statistics','-v7.3');
save('Group_Sorted_Outputs/ABCD_Group_Empirical_Statistics_Across_Models_and_Parcellations.mat','Empirical_Statistics','-v7.3');
save('Group_Sorted_Outputs/ABCD_Group_Topological_Dissimilarity.mat','Dissimilarity_Across_Models_and_Parcellations','-v7.3');
save('Group_Sorted_Outputs/ABCD_Group_Spatial_Embedding.mat','Spatial_Embedding_Across_Models_and_Parcellations','-v7.3');


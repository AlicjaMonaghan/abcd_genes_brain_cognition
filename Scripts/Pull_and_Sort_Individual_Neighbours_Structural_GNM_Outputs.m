% This script details pulling output for the neighbours structural 
% generative network model for each participant, using eta and gamma limits
% determined by the best-fitting model from a group model run across all 
% 2154 ABCD participants with a consensus-based distance-dependent network. 
% We shall sort the model energy from lowest to highest, and extract the
% best-fitting model parameters for each participant. Note that we used the
% Schaefer 100-node parcellation. Further, we shall re-run the neighbour
% structural GNM for the best-fitting parameter combination for each
% participant, in order to collect statistics throughout each run, allowing
% us to plot parameterised nodal costs and values. 
%% PART 1 - Set up the workspace.
% Load up the directory containing the individual structural GNM outputs.
clear; clc;
Structural_GNM_Outputs_RootDir = '/imaging/projects/external/abcd/analyses/Alicja/Schaefer_100_Node_Analyses/Individual_Structural_GNMs/Individual_Structural_GNM_Outputs/';
addpath(Structural_GNM_Outputs_RootDir)
cd(Structural_GNM_Outputs_RootDir);

% Adding path for the brain connectivity toolbox 
brainconnectivity_path = '/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT';
addpath(brainconnectivity_path);

%Here we are plotting 2D coordinates, spanning the eta and gamma limits,
%respectively. This grid will represent different combinations of eta and
%gamma, which shall be evaluated below (in terms of KS-statistics and
%energies).
etalimits = [-7,-.2000];
gamlimits = [-0.0667,.6000];
nruns = 75000;
[p,q] = meshgrid(linspace(etalimits(1),etalimits(2),sqrt(nruns)),...
                   linspace(gamlimits(1),gamlimits(2),sqrt(nruns)));
params = unique([p(:) q(:)],'rows');
nparams = size(params,1);

%% PART 2 - Loading Up Successfully-Processed Participants!
% Find all types of data stored in the output directory
Individual_GNM_Files = dir('Individual_Structural_neighbour_GNM_Output_*.mat');
% Now find which participants have been processed successfully
Participant_IDs = cell(size(Individual_GNM_Files));
for i = 1:length(Individual_GNM_Files)
    Participant_IDs{i} = extractBetween(Individual_GNM_Files(i).name,'Individual_Structural_neighbour_GNM_Output_',".mat");
end
% Now convert to a string for easier formatting 
Participant_IDs = string(Participant_IDs);
% And specify how many participants have been processed
nsub = length(Participant_IDs);

%% PART 2B - Update Slurm Processing to Include All Participants
% Occasionally, some participants may have been missed during processing.
% Therefore, we shall find the difference between the MATLAB data files in 
% the output directory, and those who successfully underwent structural 
% brain processing (QSIprep).

% First, find which participants have been processed by QSIprep for all 
% parcellations.
% Load up the list of successfully processed participants
cd('/imaging/projects/external/abcd/analyses/Alicja/')
Processed_Participants_Brainnetome246 = importdata('QSIprep/Successfully_Processed_Participants_brainnetome246_Parcellation.txt');
Processed_Participants_Schaefer100 = importdata('QSIprep/Successfully_Processed_Participants_schaefer100x17_Parcellation.txt');
Processed_Participants_Schaefer400 = importdata('QSIprep/Successfully_Processed_Participants_schaefer400x17_Parcellation.txt');
Processed_Participants = intersect(intersect(Processed_Participants_Schaefer100,Processed_Participants_Schaefer400),Processed_Participants_Brainnetome246);
clear Processed_Participants_Brainnetome246 Processed_Participants_Schaefer100 Processed_Participants_Schaefer400   
% Remove the headers and convert to string
Processed_Participants(1,:) = [];
Processed_Participants = string(Processed_Participants);
% And save...
writetable(table(Processed_Participants),'Structural_GNMs_Alicja_V3/Individual_GNMs/Third_Round/All_Processed_Participant_IDs.txt','Delimiter',' ','WriteVariableNames',0);
% Find the difference between the participants who have already been
% processed and those that remain to be processed.
Remaining_Participants = setdiff(Processed_Participants,Participant_IDs);
% And save
writetable(table(Remaining_Participants),'Schaefer_100_Node_Analyses/Individual_Structural_GNMs/Remaining_Participant_IDs_15th_June.txt','Delimiter',' ','WriteVariableNames',0);

%% PART 3 - Initialise Output Arrays, Pull, and Sort Output
% Also save the indices of the best simulation, so that we can re-run the
% best-performing parameter combination to retrieve parameterised
% coefficients and simulated statistics. 
ABCD_Individual_Sorted_Energy_Top_100_Simulations = zeros(nsub,100);
ABCD_Individual_Energy = zeros(nsub,nparams);
ABCD_Individual_Lowest_Energy = zeros(nsub,1);
ABCD_Individual_Best_Eta = zeros(nsub,1);
ABCD_Individual_Best_Gamma = zeros(nsub,1);
ABCD_Best_Simulation_Indices = zeros(nsub,1);

% Now pull and sort output for all participants.
for sub = 1:nsub
    try 
        fprintf('Pulling output for %s! \n', Participant_IDs(sub));
        % Start the timer...
        tic 
        % Load the output file for each participant
        Output_File_Name = sprintf('Individual_Structural_neighbour_GNM_Output_%s.mat',Participant_IDs(sub));
        Loaded_Output = load(Output_File_Name);
        Loaded_Output = Loaded_Output.Individual_Structural_Neighbours_Output;
        % Load the energy for the model
        ABCD_Individual_Energy(nsub,:) = Loaded_Output.Energy;
        % Find the lowest 10% energy for the participant
        [~,Top_100_Runs_Indices] = sort(Loaded_Output.Energy);
        Top_100_Runs_Indices = Top_100_Runs_Indices(1:100);
        ABCD_Individual_Sorted_Energy_Top_100_Simulations(sub,:) = Loaded_Output.Energy(Top_100_Runs_Indices);
        % Now find the lowest energy, and the associated eta and gamma
        % parameters!
        [~,minimum_idx] = min(Loaded_Output.Energy);
        ABCD_Individual_Lowest_Energy(sub,1) = Loaded_Output.Energy(minimum_idx);
        % Extracting the associated eta
        ABCD_Individual_Best_Eta(sub,1) = params(minimum_idx,1);
        % And the associated gamma
        ABCD_Individual_Best_Gamma(sub,1) = params(minimum_idx,2);
        % Save the minimum index for each participant
        ABCD_Best_Simulation_Indices(sub,1) = minimum_idx;
        % Clear unnecessary variables to save space
        clear Loaded_Output 
        % Stop the timer and update the user
        Time_Elapsed = toc;
        fprintf('Successfully pulled output for %s in %d seconds! \n', Participant_IDs(sub), round(Time_Elapsed,2));
    catch 
    end
end

% Save these outputs into a new sorted output directory
save('Sorted_Outputs/ABCD_Individual_Sorted_Energy_Top_100_Simulations.mat','ABCD_Individual_Sorted_Energy_Top_100_Simulations','-v7.3');    
save('Sorted_Outputs/ABCD_Individual_Lowest_Energy.mat', 'ABCD_Individual_Lowest_Energy', '-v7.3');
save('Sorted_Outputs/ABCD_Individual_Best_Eta.mat','ABCD_Individual_Best_Eta','-v7.3');
save('Sorted_Outputs/ABCD_Individual_Best_Gamma.mat','ABCD_Individual_Best_Gamma','-v7.3');
save('Sorted_Outputs/ABCD_Best_Simulation_Indices.mat','ABCD_Best_Simulation_Indices','-v7.3');
% And save the participant IDs of those processed
writetable(table(Participant_IDs), 'Sorted_Outputs/Processed_Participant_IDs_Accompanying_Sorted_Outputs.txt','Delimiter',' ','WriteVariableNames',0);

% Now clear the majority of the above variables to save space
clear ABCD_Individual_Sorted_Energy_Top_100_Simulations p q 

%% PART 3B - Visualise Energy Landscape and Optimal Parameters %%
Landscape = squeeze(mean(reshape(ABCD_Individual_Energy, [2153 273 273]),1));
imagesc(Landscape)
xlabel('\fontsize{30} \eta','Interpreter','tex', 'FontWeight','bold', 'Position',[150 270])
ylabel('\fontsize{30} \gamma','Interpreter','tex', 'FontWeight','bold','Position',[-30 150], 'Rotation', 360)
xticks([2 270]); 
xticklabels({round(min(ABCD_Individual_Best_Eta(:)),3), ...
    round(max(ABCD_Individual_Best_Eta(:)),3)});
set(gca, 'TickLength', [0 0])
a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'FontSize', 15)
yticks([2 270]); yticklabels({round(min(ABCD_Individual_Best_Gamma(:)),3), ...
    round(max(ABCD_Individual_Best_Gamma(:)),3)});
axis square
box off
hold on
c = colorbar;
set(c, 'TickLength', 0, 'FontSize', 15)
set(get(c, 'ylabel'), 'string', 'Energy', 'FontSize', 20, 'FontWeight', 'bold', 'Rotation', 90)
% Rescale the best individual eta and gamma values to the limits of the
% energy landscape plot. The x and y limits are the same.
x1 = xlim;
ABCD_Individual_Best_Parameters_Scaled = horzcat(rescale(ABCD_Individual_Best_Eta, x1(1), x1(2)), rescale(ABCD_Individual_Best_Gamma, x1(1), x1(2)));
scatter(ABCD_Individual_Best_Parameters_Scaled(:,1), ABCD_Individual_Best_Parameters_Scaled(:,2), ...
    80, 'x', 'MarkerEdgeColor', 'black', 'LineWidth', 1.5)
hold off
saveas(gcf, ['/imaging/projects/external/abcd/analyses/Alicja/Schaefer_100_' ...
    'Node_Analyses/Individual_Structural_GNMs/Visualisation/Narrow_' ...
    'Landscape_Optimal_Parameters.jpg'], 'jpg')

%% PART 3B - For the Best-Performing Individual Parameter Combination, Calculate 
%  Simulated Statistics, Alongside Node-Wise and Edge-Wise Costs and Values 
%  Within Each Model %%

% To calculate spatial embedding, we need to obtain simulated
% statistics for the best-performing parameter combination for each
% participant. Therefore, for each participant, using the indices of their
% best-performing parameter combination, we shall re-run the neighbour
% generative model, and save the outputs. Further, for the best-performing 
% pair of parameters for each participant, we shall re-run the neighbour 
% structural model, and collect node-wise and edge-wise statistics every 
% time a new edge is added. Since this is time-consuming and we have many
% participants, we created a function called 'Run_Nodal_Parameterised_and_
% Simulated_Statistics', output of which we shall collect and sort here.

% Find the participants whose statistics have been collected.
Individual_Nodal_Parameterised_and_Nodal_Statistics = dir('/imaging/projects/external/abcd/analyses/Alicja/Schaefer_100_Node_Analyses/Individual_Structural_GNMs/Nodal_Parameterised_and_Simulated_Statistics/sub-*.mat');
% Extract the participant IDs
Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs = cell(size(Individual_Nodal_Parameterised_and_Nodal_Statistics));
for sub = 1:length(Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs)
    Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs{sub} = Individual_Nodal_Parameterised_and_Nodal_Statistics(sub).name(1:end-4);
end
Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs = string(Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs);
% And save the above!
writetable(table(Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs), '/imaging/projects/external/abcd/analyses/Alicja/Schaefer_100_Node_Analyses/Individual_Structural_GNMs/Nodal_Parameterised_and_Simulated_Statistics/Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs.txt','Delimiter',' ','WriteVariableNames',0);

% Initialise output arrays:
% Collect all simulated statistics for all participants.
All_Simulated_Statistics_Optimal_Parameters = cell(10,length(Individual_Nodal_Parameterised_and_Nodal_Statistics));
% Node-wise parameterised cost (D)
Neighbour_Parameterised_Nodal_Cost = zeros(length(Individual_Nodal_Parameterised_and_Nodal_Statistics),100,100);
% Node-wise parameterised value (K)
Neighbour_Parameterised_Nodal_Value = zeros(length(Individual_Nodal_Parameterised_and_Nodal_Statistics),100,100);

% Now loop over each participant, and extract their statistics, appending
% it to the arrays above.
for sub = 1:length(Individual_Nodal_Parameterised_and_Nodal_Statistics)
    % Load up the output for this participant
    Loaded_Output = load(sprintf(['/imaging/projects/external/abcd/analyses/Alicja/' ...
        'Schaefer_100_Node_Analyses/Individual_Structural_GNMs/Nodal_Parameterised_' ...
        'and_Simulated_Statistics/%s'],Individual_Nodal_Parameterised_and_Nodal_Statistics(sub).name));
    Loaded_Output = Loaded_Output.Nodal_Parameterised_and_Simulated_Statistics;
    % Assign the simulated statistics to the output
    for i = 1:10
        All_Simulated_Statistics_Optimal_Parameters{i,sub} = Loaded_Output.Simulated_Statistics{i,1};
    end
    % Assign the node-wise parameterised nodal cost
    Neighbour_Parameterised_Nodal_Cost(sub,:,:) = squeeze(Loaded_Output.Parameterised_Nodal_Cost);
    % Assign the node-wise parameterised nodal value
    Neighbour_Parameterised_Nodal_Value(sub,:,:) = squeeze(Loaded_Output.Parameterised_Nodal_Value);
end

% We now need to calculate the mean neighbourhood parameterised nodal costs
% and value across participants!
Mean_Neighbour_Parameterised_Nodal_Cost = sum(squeeze(mean(Neighbour_Parameterised_Nodal_Cost)),2);
Mean_Neighbour_Parameterised_Nodal_Value = sum(squeeze(mean(Neighbour_Parameterised_Nodal_Value)),2);
% Save the outputs!
cd('Schaefer_100_Node_Analyses/Individual_Structural_GNMs/')
save('Sorted_Outputs/All_Simulated_Statistics_Optimal_Parameters.mat','All_Simulated_Statistics_Optimal_Parameters','-v7.3')
save('Sorted_Outputs/Neighbour_Parameterised_Nodal_Cost.mat','Neighbour_Parameterised_Nodal_Cost','-v7.3');
save('Sorted_Outputs/Neighbour_Parameterised_Nodal_Value.mat','Neighbour_Parameterised_Nodal_Value','-v7.3');
save('Sorted_Outputs/Mean_Neighbour_Parameterised_Nodal_Value_Across_Participants.mat','Mean_Neighbour_Parameterised_Nodal_Value','-v7.3');
save('Sorted_Outputs/Mean_Neighbour_Parameterised_Nodal_Cost_Across_Participants.mat','Mean_Neighbour_Parameterised_Nodal_Cost','-v7.3');

%% PART 4 - Collect Output from Each Participant's Best Simulation!
% To find how similar the simulated and observed connectomes are, we shall
% compute the spatial embedding for each, as well as correlate the
% different network statistics. Therefore, for each participant's 
% best-performing simulation, we need to grab the simulated and observed statistics, 

% Now load the coordinates for the Schaefer 400-node parcellation!
Coordinates = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x17_1mm_info.mat');
D = squareform(pdist([Coordinates.schaefer100x17_1mm_info.x_mni Coordinates.schaefer100x17_1mm_info.y_mni Coordinates.schaefer100x17_1mm_info.z_mni]));
nroi = length(D);

% Load up each participant's target connectome.
All_Participant_Targets = load('/imaging/projects/external/abcd/analyses/Alicja/Structural_GNMs_Alicja_V3/Thresholding_and_Binarisation/ABCD_Individual_Target_Connectomes_schaefer100x17_Parcellation.mat');
All_Participant_Targets = All_Participant_Targets.ABCD_Thresholded_27_Streamlines;
% And load the participant IDs which accompany this!
Passed_All_Parcellations = importdata('/imaging/projects/external/abcd/analyses/Alicja/Structural_GNMs_Alicja_V3/Passed_All_Parcellations.txt');
% Remove the participant whose processing failed
[~,idx] = intersect(Passed_All_Parcellations, Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs);
All_Participant_Targets = All_Participant_Targets(idx,:,:);
% Finally, load up the seed common to all participants
Seed = load('/imaging/projects/external/abcd/analyses/Alicja/Structural_GNMs_Alicja_V3/Thresholding_and_Binarisation/ABCD_Seed_Network_schaefer100x17_Parcellation.mat');
Seed = Seed.ABCD_Seed_Network;

% Finally, initialise an output variable to hold each participant's
% spatial embedding values, alongside simulated and observed statistics.
Spatial_Embedding_All_Individual_Neighbours_Models = zeros(length(Individual_Nodal_Parameterised_and_Nodal_Statistics),1);
All_Observed_Statistics = cell(length(Individual_Nodal_Parameterised_and_Nodal_Statistics),10);
Correlation_Observed_Simulated_Statistics = zeros(length(Individual_Nodal_Parameterised_and_Nodal_Statistics),7);
% For the below, the first entry of the second dimension is the observed
% non-zero edge length, whilst the second part of the second dimension is
% the simulated non-zero edge length.
Non_Zero_Edge_Length = cell(length(Individual_Nodal_Parameterised_and_Nodal_Statistics),2);

% Now, run the generative models, and extract outputs!
for sub = 1:length(Individual_Nodal_Parameterised_and_Nodal_Statistics)
    % Load the observed statistics for this participant.
    Output_Observed_Statistics = load(sprintf('Individual_Structural_neighbour_GNM_Observed_Statistics_%s.mat',Participant_IDs(sub)));
    Observed_Statistics = Output_Observed_Statistics.x;
    % Assign to the output variable!
    for i = 1:10
        All_Observed_Statistics{sub,i} = Observed_Statistics{i,1};
    end
    % Now load the individual target connectome!
    Target = squeeze(All_Participant_Targets(sub,:,:));
    % Now calculate the total distance from each node, after re-calculating
    % the total edge length of each node.
    All_Observed_Statistics{sub,4} = sum(D.*(Target>0))';
    Node_Edge_Length_Observed = nroi;
    for n = 1:nroi
        s = All_Observed_Statistics{sub,4}(n,:);
        Node_Edge_Length_Observed(n) = sum(norm(n-s));
    end
    % And assign to the correct part of the output data frame
    All_Observed_Statistics{sub,4} = Node_Edge_Length_Observed;
    % Finally, collect the length of non-zero edges
    Non_Zero_Edge_Length{sub,1} = D(triu(Target,1) > 0);
    % Clear variables we don't require!
    clear Output_Observed_Statistics Observed_Statistics Node_Edge_Length_Observed

    % Re-run the best-performing simulation for this participant, in order
    % to calculate the total simulated edge length.
    n = length(D);
    m = nnz(Target)/2;
    b = zeros(m,nparams);
    modelvar = [{'powerlaw'},{'powerlaw'}]; 
    epsilon = 1e-5;
    Individual_Parameters = params(ABCD_Best_Simulation_Indices(sub,1),:);
    B = generative_model(Seed,D,m,"neighbors",modelvar,Individual_Parameters,epsilon);
    %Finding the size of the second dimension of the generative model 
    nB = size(B,2); 
 for iB = 1:nB %For the number of columns/parameters in the generative network model
     b = zeros(n); %Create a matrix of zeros for the cardinality of the network (size of Connectivity Dimensions)
     b(B(:,iB)) = 1; %For iB indexing the parameter of the generative network model, find the associated column in the cardinality network
     b = b + b'; %Represent b as an adjacency matrix 

     % Calculate the synthetic model properties!
     All_Simulated_Statistics_Optimal_Parameters{4,sub}  = sum(D.*(b>0))';
     Non_Zero_Edge_Length{sub,2} = D(triu(b,1) > 0); 
 end
    % For the simulated statistics, also calculate the simulated node edge
    % length!
    Node_Edge_Length_Simulated = nroi;
        for n = 1:nroi
            s = All_Simulated_Statistics_Optimal_Parameters{4,sub}(n,:);
            Node_Edge_Length_Simulated(n) = sum(norm(n-s));
        end
    % And assign to the correct part of the output data frame.
    All_Simulated_Statistics_Optimal_Parameters{4,sub} = Node_Edge_Length_Simulated;
    % Clear variables
    clear Node_Edge_Length_Simulated

    % And now correlate only the local statistics to retrieve the spatial
    % embedding!
    for i = 1:6
        Spatial_Embedding_Correlation = corrcoef(All_Observed_Statistics{sub,i},All_Simulated_Statistics_Optimal_Parameters{sub,i});
        Spatial_Embedding_All_Individual_Neighbours_Models(sub,:) = Spatial_Embedding_Correlation(1,2);
    end
    % Clear unnecessary variables
    clear All_Measures_Correlation Spatial_Embedding_Correlation
        
end

% Append the non-zero edge length to each of the relevant data frames!
All_Observed_Statistics = horzcat(Non_Zero_Edge_Length(:,1), All_Observed_Statistics);
All_Simulated_Statistics_Optimal_Parameters = horzcat(All_Simulated_Statistics_Optimal_Parameters',Non_Zero_Edge_Length(:,2));

% Now save the synthetic matrices (PART 4) and the spatial embedding values
% (PART 5).
cd('Schaefer_100_Node_Analyses/Individual_Structural_GNMs/')
save('Individual_Neighbours_Simulated_Statistics_Best_Simulation.mat','All_Simulated_Statistics_Optimal_Parameters');
save('Individual_Neighbours_Observed_Statistics_Best_Simulation.mat','All_Observed_Statistics');
save('Individual_Neighbours_Models_Spatial_Embedding.mat', 'Spatial_Embedding_All_Individual_Neighbours_Models');

%% PART 6 - Correlate Individual Observed and Simulated Statistics! %%
% Calculate the within-participant mean for each measure
Measures_Indices = [1,2,3,4,5,6,8,10,11];
Within_Participant_Averages_Observed_Statistics = zeros(nsub,9);
Within_Participant_Averages_Simulated_Statistics = zeros(nsub,9);
for sub = 1:nsub
    for Measure = Measures_Indices
        Within_Participant_Averages_Observed_Statistics(sub,Measure) = mean(All_Observed_Statistics{sub,Measure});
        Within_Participant_Averages_Simulated_Statistics(sub,Measure) = mean(All_Simulated_Statistics_Optimal_Parameters{sub,Measure});
    end
end

% Calculate across-participant means for each measure
Across_Participants_Averages_Observed_Statistics = zeros(nroi,9);
Across_Participants_Averages_Simulated_Statistics = zeros(nroi,9);

for Measure = Measures_Indices
    Concatenated_Measure_Observed = squeeze(cat(3,All_Observed_Statistics{:,Measure}));
    Across_Participants_Averages_Observed_Statistics(:,Measure) = mean(Concatenated_Measure_Observed,2);
    Concatenated_Measure_Simulated = squeeze(cat(3,All_Simulated_Statistics_Optimal_Parameters{:,Measure}));
    Across_Participants_Averages_Simulated_Statistics(:,Measure) = mean(Concatenated_Measure_Simulated,2);
end

% And correlate each observed and simulated measure!
% Within_Participant_Correlations(:,1) is the correlation coefficient,
% whilst Within_Participant_Correlations(:,2) is the p-value for the null
% hypothesis of no significant association.
Within_Participant_Correlations = zeros(8,2);
for Measure = Measures_Indices
    [R,P] = corrcoef(Across_Participants_Averages_Observed_Statistics(:,Measure), Across_Participants_Averages_Simulated_Statistics(:,Measure));
    Within_Participant_Correlations(Measure,1) = R(1,2);
    Within_Participant_Correlations(Measure,2) = P(1,2);
end

% Save the above outputs!
save('Sorted_Outputs/Within_Participant_Averaged_Observed_Statistics.mat','Within_Participant_Averages_Observed_Statistics','-v7.3');
save('Sorted_Outputs/Within_Participant_Averaged_Simulated_Statistics.mat','Within_Participant_Averages_Simulated_Statistics','-v7.3');
save('Sorted_Outputs/Across_Participants_Averages_Observed_Statistics.mat','Across_Participants_Averages_Observed_Statistics','-v7.3');
save('Sorted_Outputs/Across_Participants_Averages_Simulated_Statistics.mat','Across_Participants_Averages_Simulated_Statistics','-v7.3');

%% PART 7 - Find Similarity Between Synthetic and Observed Connectomes.
% As described above, we shall compare the synthetic and observed
% connectomes using a support vector machine classifier.

% First, concatenate the simulated and observed connectomes, after
% reshaping each into row vectors. 
Simulated = reshape(Asynthall,[1127 160000]);
Observed = reshape(Target_Connectomes,[1127 160000]);
% The below shows the simulated and observed connectomes appened
% end-to-end, where the first half of the first dimension corresponds to
% individual synthetic networks, and the second half of the first dimension
% corresponds to individual observed networks.
Simulated_and_Observed_Connectomes = cat(1,Simulated,Observed);
% Now create the Simulated and Observed classes!
Simulated_Class_Label = repmat('Simulated',1127,1);
Observed_Class_Label = repmat('Observed',1127,1);
Class_Labels = vertcat(string(repmat('Simulated',1127,1)), string(repmat('Observed',1127,1)));
% And fit the SVM!
Simulated_Observed_SVM = fitcsvm(Simulated_and_Observed_Connectomes, Class_Labels);



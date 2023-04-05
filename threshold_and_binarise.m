% This script details structural connectome binarisation, thresholding,  
% setting the seed and target networks for generative network modelling 
% (GNM), for a stratified subset of 2154 participants from the ABCD 
% Baseline Release. To remove spurious connections, we first conduct 60%
% consensus thresholding. To create the group GNM target, we conduct
% distance-dependent consensus thresholding (Betzel et al., 2019). To
% create the individual-level GNM target, we threshold individual
% connectomes at 27 streamlines. To create the seed networks for group- and
% individual-level analyses, we threshold individual connectomes at a 95%
% consensus. Note that this code is for the Schaefer 100-node (17-network)
% parcellation, for which most of the analysis uses. Correspondence should 
% be directed to Alicja Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk

% STEPS:
% 1. Preparing the work space and loading data.
% 2. 60% consensus thresholding across all participants.
% 3. Distance-dependent consensus group connectome: Group-level GNM target.
% 4. Thresholding individual connectomes at 27 streamlines:
% Individual-level GNM target.
% 5. Seed network (for all group- and individual-level GNMs) - 
% Supplementary Figure 2.

%% PART 1 - PREPARE THE WORKSPACE AND LOAD DATA %%
clear; clc;
% Set working directory to toolboxes, code, and data. <-- SET THIS TO WHERE
% YOU HAVE SAVED THIS DIRECTORY.
cd('abcd_genomic_variation_structural_generative_mechanisms_open/');
% Add path to the brain connectivity toolbox (Rubinov and Sporns, 2010),
% and the distance-dependent consensus thresholding procedure (Betzel and
% colleagues, 2019)
addpath('toolboxes_and_functions/2019_03_03_BCT');
addpath('toolboxes_and_functions/distanceDependent');
% Load the QSIprep structural connectomes for the 2154 participants with
% reconstructed connectomes across all 3 parcellations (Schaefer 100-node
% 17-network, Brainnetome 246-node, and Schaefer 400-node 17-network).
structural_connectomes = load('data/qsiprep_structural_connectomes.mat');
structural_connectomes = structural_connectomes.structural_connectomes;
% Load the meta-data for the Schaefer 100-node 17-network parcellation.
schaefer100_metadata = load('data/schaefer100x17_1mm_info.mat');
% Load the coordinates for the parcellation
schaefer100_coordinates = [schaefer100_metadata.schaefer100x17_1mm_info.x_mni,...
    schaefer100_metadata.schaefer100x17_1mm_info.y_mni,... 
    schaefer100_metadata.schaefer100x17_1mm_info.z_mni];
% Calculate the euclidean distance between the coordinates.
D = squareform(pdist(schaefer100_coordinates));
% Set the number of regions for this parcellation.
nroi = 100;
% And set the number of participants.
nsub = length(structural_connectomes);

%% PART 2 - 60% CONSENSUS THRESHOLDING ACROSS ALL PARTICIPANTS %%
% To remove spurious connections, and in line with previous connectomic 
% studies, we threshold all individual connectomes to include only 
% connections common to at least 60% of participants.

% First, remove self-connections.
for sub = 1:nsub
    structural_connectomes(sub,:,:) = squeeze(structural_connectomes(sub,:,:)) - diag(diag(squeeze(structural_connectomes(sub,:,:))));
end

% Set the consensus threshold
Consensus_Threshold = .60;
t = floor(nsub*Consensus_Threshold);
% Find non-zero elements
k = structural_connectomes~=0;
u = squeeze(sum(k,1));
% Identify indices to remove
ind = u<=t;

% Initialise matrices for the 60% consensus thresholded matrices.
abcd_60_percent_consensus_thresholded_structural_connectomes = zeros(nsub, nroi, nroi);
% Apply the consensus thresholding to each participant.
for sub = 1:nsub
    A = squeeze(structural_connectomes(sub,:,:));
    % Remove edges
    A(ind) = 0;
    % Assign to output array
    abcd_60_percent_consensus_thresholded_structural_connectomes(sub,:,:) = A;
end

%% PART 3 - CREATE A DISTANCE-DEPENDENT CONSENSUS GROUP CONNECTOME %%
% This will act as the target network for the group model. We shall be using
% a distance-dependent consensus-group algorithm developed by Betzel and
% colleagues (2019) demonstrated to preserve the distribution of shorter
% edges. The input to this consensus network will be the weighted
% individual connectomes thresholded at a 60% consensus.

% Set the hemisphere IDs i.e. 1 signifying the left hemisphere, and 2
% signifying the right hemisphere.
hemiid = zeros(nroi,1);
hemiid(1:nroi/2,1) = 1;
hemiid(nroi/2:end,1) = 2;

% Transpose the connectome so that dimensions align i.e. nroi x nroi x nsub
abcd_60_percent_consensus_thresholded_structural_connectomes = permute(abcd_60_percent_consensus_thresholded_structural_connectomes,[2 3 1]);
% Calculate the group distance-dependent consensus-thresholded connectome
[group_distance_binarised_structural_connectome, ~] = fcn_group_bins(abcd_60_percent_consensus_thresholded_structural_connectomes,D,hemiid,100);

%% PART 4 - BINARISE INDIVIDUAL CONNECTOMES %%
% We inputted the weighted individual connectomes into the group consensus,
% where each individual connectome only included connections present in at
% least 60% of participants. However, to create our seed and targets for
% the individual models, we need binarised networks. Therefore, we shall
% first threshold the individual connectomes at 27 streamlines to produce
% the individual target networks.

% Now initialise an empty dataframe for the binarised connectomes for all
% participants. Also initialise an output array for the 4 graph theory
% metrics we'll visualise across participants in Figure 1d-g i.e. degree,
% edge length, betweenness-centrality, and clustering coefficient.
abcd_thresholded_27_streamlines = zeros(nsub, nroi, nroi);
abcd_thresholded_27_streamlines_density = zeros(nsub,1);
abcd_individual_connectome_properties = zeros(nsub,nroi,4);
threshold = 27;
% Loop through participants. 
for sub = 1:nsub
    % Extract the individual's connectome weighted connectome thresholded
    % at 60% consensus.
    W = squeeze(abcd_60_percent_consensus_thresholded_structural_connectomes(:,:,sub)); 
    % Calculate the absolute threshold
    W_thr = threshold_absolute(W,threshold);  
    % Calculate the density and append to output array
    abcd_thresholded_27_streamlines_density(sub,1) = density_und(W_thr)*100;  
    B = zeros(size(W_thr));
    % Binarise
    B(find(W_thr)) = 1;  
    % Append to output array
    abcd_thresholded_27_streamlines(sub,:,:) = B;
    % For the thresholded connectome, calculate degree...
    abcd_individual_connectome_properties(sub,:,1) = sum(B,2);
    % And edge length...
    abcd_individual_connectome_properties(sub,:,2) = sum(D.*(B>0));
    % And betweenness-centrality...
    abcd_individual_connectome_properties(sub,:,3) = betweenness_bin(B);
    % And clustering coefficient...
    abcd_individual_connectome_properties(sub,:,4) = clustering_coef_bu(B);
end
% Update user about the mean density across participants.
fprintf(['For the Schaefer 100-node parcellation, the mean density across participants is %.2f percent, ' ...
    'with standard deviation %.2f percent.\n'], mean(abcd_thresholded_27_streamlines_density,1), ...
    std(abcd_thresholded_27_streamlines_density,1));

%% PART 5 - COMPUTE AND VISUALISE THE SEED NETWORK %%
% For both the group and individual models, the seed shall be the
% connections present across at least 95% of participants. 
seed_network = zeros(nroi,nroi);
connections = squeeze(mean(abcd_thresholded_27_streamlines,1));
index = find(connections > .95);
seed_network(index) = 1;
seed_network = upper(seed_network);
% Update user about the density
fprintf('The seed network has a density of %.2f percent.\n', density_und(seed_network)*100);
% Save!
save('data/example_seed_gnm.mat','seed_network');
% Now visualise the seed network - this will be SUPPLEMENTARY FIGURE 2.
imagesc(seed_network);
title('Seed Network', 'FontSize',12, 'FontWeight','bold');
subtitle('157 common connections across at least 95% of participants', 'FontSize',12);
xlabel('Node','FontSize',10);
ylabel('Node', 'FontSize',10);
ax = gca;
set(gca, 'XTick', [], 'YTick', []);
% Add colourbar
c = colorbar('southoutside', 'Ticks',[0,0.5,1], ...
    'TickLabels',{'Absent', 'Common Connections', 'Present'}, 'TickLength', [],'FontSize', 10);
exportgraphics(gca,'seed_network.png','Resolution',600)
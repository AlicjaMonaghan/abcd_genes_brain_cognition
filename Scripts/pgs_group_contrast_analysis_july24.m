% This script compares the topological and graph theory distributions for
% the top and bottom 10% polygenic score (PGS) groups. Written by Danyal
% Akarca and Alicja Monaghan, University of Cambridge.
%% add pre-requisites
clear; clc;
% change directory
cd('/imaging/astle/users/da04/Postdoc/pgs_dissim/outputs');
% load generative models
load('group_1_generative_model.mat'); group_1 = output.networks; group_1_prob = output.probabilities;
load('group_2_generative_model.mat'); group_2 = output.networks; group_2_prob = output.probabilities;
% load seed
load('/imaging/astle/users/da04/Postdoc/pgs_dissim/seed.mat');
% load coordinates
load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x17_1mm_info');
coordinates = [schaefer100x17_1mm_info.x_mni schaefer100x17_1mm_info.y_mni schaefer100x17_1mm_info.z_mni]; 
euclidean = squareform(pdist(coordinates));
% form statistics labels
stat_labels = string({'Clustering','Betweenness','Edge length','Efficiency','Modularity','Communicability'});
% for dissimilarity labels
dissim_labels = string({'\SigmaEuclidean','Consistency'});
% addpath bct
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath cohens d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');

%% compute node-wise statistics for each network
% set nnodes
nnode = 100;
% set number of runs
nrun = 1000;
% set number of parameters
nparams = 2;
% set number of global statistics measured
ngmeasures = 6;
% set output
output = [];
output(:,:,1) = group_1;
output(:,:,2) = group_2;
% set euclidean
euclidean = euclidean;
% compute node wise statistics
global_statistics = zeros(nrun,nparams,ngmeasures);
% loop over repetitions
for rep = 1:nrun;
    % loop over parameter combinations
    for parametercomb = 1:nparams;
        % form network
        A = zeros(nnode,nnode);
        ind = squeeze(output(rep,:,parametercomb));
        A(ind) = 1;
        A = A + A';
        % compute node statisitcs
        % mean clustering coefficient
        group_statistics(rep,parametercomb,1) = mean(clustering_coef_bu(A));
        % mean betweenness centrality
        group_statistics(rep,parametercomb,2) = mean(betweenness_bin(A));
        % total edge lengths
        group_statistics(rep,parametercomb,3) = sum(euclidean.*A,'all');
        % global efficiency
        group_statistics(rep,parametercomb,4) = efficiency_bin(A);
        % modularity q
        [~,group_statistics(rep,parametercomb,5)] = modularity_und(A);
        % communicability
        group_statistics(rep,parametercomb,6) = mean(expm(A),'all');
        % display
        fprintf('parametercomb %g rep %g node statistics computed\n',parametercomb,rep);
    end
end
% save('group_statistics.mat');
%% compute dissimilarity within parameter combination
% load data
load('group_statistics.mat');
% assign
globalstatistics = group_statistics;
% set nnodes
nnode = 100;
% set number of runs
nrun = 1000;
% set number of parameters
nparams = 2;
% set number of connections
nconnect = 331;
% set number of dissimilairty measures computed
ndissim = 2;
% initialise
dissimilarity = zeros(nparams,ndissim);
hdistances = zeros(nparams,nrun,nrun);
cdistances = zeros(nparams,nrun,nrun);
% all combinations
com = combnk(1:nrun,2);
% compute dissimilarity between same parameter runs
for parameterscomb = 1:nparams;
    % get all the runs for this parameter
    data = squeeze(globalstatistics(:,parameterscomb,:)); % check which dimension is parameter, which is repetition
    % normalize the data
    data = normalize(data);
    % calculate the high-dimensional distance
    hdist = squareform(pdist(data));
    % keep
    hdistances(parameterscomb,:,:) = hdist;
    % measure of similarity 1: sum euclidean
    dissimilarity(parameterscomb,1) = sum(hdist,'all');
    % sum up the shared connectivity across runs
    connections = squeeze(cdistances(parameterscomb,:,:));
    mask = triu(true(size(connections)),1);
    dissimilarity(parameterscomb,2) = mean(connections(mask));
    % display
    disp(sprintf('parametercomb %g dissimilarity computed',parameterscomb));
end
embedding_topological_dissimilarity = struct;
embedding_topological_dissimilarity.embedding_matrices = cdistances;
embedding_topological_dissimilarity.topological_matrices = hdistances;
embedding_topolgoical_dissimilarity.dissimilarity = dissimilarity;
% save('embedding_topological_dissimilarity.mat');
%% look at the difference between the matrices and probabilities
load('embedding_topological_dissimilarity.mat');
% embedding
U = triu(ones(1000),1);
group_1_embedding = squeeze(embedding_topological_dissimilarity.embedding_matrices(1,:,:)); group_1_embedding = group_1_embedding(find(U));
group_2_embedding = squeeze(embedding_topological_dissimilarity.embedding_matrices(2,:,:)); group_2_embedding = group_2_embedding(find(U));
[h,p] = ttest(group_1_embedding,group_2_embedding);
d = computeCohen_d(group_1_embedding,group_2_embedding);
fprintf('embedding dissimilarity p=%.3g, d=%.3g\n',p,d);
% topological
U = triu(ones(1000),1);
group_1_topological = squeeze(embedding_topological_dissimilarity.topological_matrices(1,:,:)); group_1_topological = group_1_topological(find(U));
group_2_topological = squeeze(embedding_topological_dissimilarity.topological_matrices(2,:,:)); group_2_topological = group_2_topological(find(U));
[h,p,~,stats] = ttest(group_1_topological,group_2_topological);
d = computeCohen_d(group_1_topological,group_2_topological);
fprintf('topological dissimilarity p=%.3g, d=%.3g\n',p,d);
% probabilities
group_1_cv = std(group_1_prob,[],2)./mean(group_1_prob,2);
group_2_cv = std(group_2_prob,[],2)./mean(group_2_prob,2);
[h,p] = ttest(group_1_cv,group_2_cv);
d = computeCohen_d(group_1_cv,group_2_cv);
fprintf('coefficienct of variation p_ij p=%.3g, d=%.3g\n',p,d);
% Compare the graph theory metrics included in the topological
% dissimilarity calculation for the low and high PGS groups.
for stat_idx = 1:length(stat_labels)
    group_1 = squeeze(group_statistics(:, 1, stat_idx));
    group_2 = squeeze(group_statistics(:, 2, stat_idx));
    d = computeCohen_d(group_1, group_2);
    [~, p] = ttest(group_1, group_2);
    fprintf('%s: p=%.3g, d=%.3g\n', stat_labels{stat_idx}, p, d);
end

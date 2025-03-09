% This code compares the developmental trajectory of global efficiency for
% participants belonging to either the top or bottom 10% of the PGS
% distribution. To probe the effect of changing parameters on stochasticity
% and network architecture, we conduct a simulation experiment where we
% simulate connectivity at decreasing deciles of the optimal group-level
% eta and gamma estimates. These analyses are inspired by Sofia Carozza's
% 2023 Developmental Psychobiology paper "Early adversity changes the
% economic conditions of mouse structural brain organization". Written by
% Alicja Monaghan in February 2025.

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
%% Part 2 - Generative network models
% The first two eta-gamma parameter combinations will be for the bottom and
% top 10% PGS distribution, for which we'll collect step-wise global
% efficiency estimates. The remainder will be the best-fitting group-level
% set of parameters in deciles, for which we'll collect global efficiency
% for the final simulated networks.
eta_params = [-3.0610, -2.9630];
gamma_params = [.2138, .2144];
nparams = length(eta_params);
nruns = 1000;
n = length(seed);
% Find the number of steps the simulations will take
nsteps = m-mseed;
% Initialise to hold the step-wise global efficiency for the PGS contrasts
neighbors_global_efficiency = cell(nruns, nparams);
% Rescale the number of simulations, including the seed network, to between
% 0 and 1.
rescaled_steps = rescale(1:(nsteps+1));
%% Part 2A - Collect step-wise global efficiency for PGS contrasts
for iparams = 1:nparams
    for run = 1:nruns
        Kseed = (seed*seed).*~eye(n);
        fprintf('running run %g of %g for parameter combination %g...\n',run,nruns, iparams);
        stepwise_networks = fcn_neighbors(seed,Kseed,D,m,eta_params(iparams),gamma_params(iparams),epsilon);
        % Initialise array to hold step-wise global efficiency for this run
        % and parameter combination. Add 1 to nsteps to account for the
        % seed network.
        stepwise_global_efficiency = zeros(nsteps+1,1);
        % Calculate global efficiency at each step of the simulation
        for step=1:(nsteps+1)
            stepwise_global_efficiency(step,1) = efficiency_bin(squeeze(stepwise_networks(step,:,:)));
        end
        % Assign back to cell!
        neighbors_global_efficiency{run, iparams} = stepwise_global_efficiency;
    end
    % Initialise an array to hold the global efficiency values across all
    % steps and runs, for this parameter combination.
    global_efficiency_array = zeros(nsteps+1, nruns);
    for run=1:nruns
        global_efficiency_array(:, run) = cell2mat(neighbors_global_efficiency(run, iparams));
    end
    % Initialise an array to hold the mean and standard deviation of global
    % efficiency at each step
    global_efficiency_summary_array = zeros(nsteps+1, 2);
    % For each step, calculate the mean and standard deviation
    for step = 1:(nsteps+1)
        global_efficiency_summary_array(step, 1) = mean(global_efficiency_array(step,:));
        % If we're examining the first step, this is the seed network, for
        % which the standard deviation is 0
        if step == 1
            global_efficiency_summary_array(step, 2) = 0;
        else
            global_efficiency_summary_array(step, 2) = std(global_efficiency_array(step,:));
        end
    end
    % Create a summary table
    global_efficiency_group_table = array2table(global_efficiency_summary_array,...
        'VariableNames',{'mean', 'std'});
    % Add the scaled steps
    global_efficiency_group_table.("rescaled_steps") = transpose(rescaled_steps);
    % And save
    if iparams == 1
        filename = 'data/stepwise_global_efficiency_table_low_PGS.csv';
    else
        filename = 'data/stepwise_global_efficiency_table_high_PGS.csv';
    end
    writetable(global_efficiency_group_table, filename, 'Delimiter',',');
end
%% Part 2B - Collect group global efficiency of final simulated networks
deciles = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
synthetic_connectome_global_efficiency = zeros(length(deciles), nruns);
% We're using the optimal eta and gamma parameters at the group level, with
% the same seed and target as before.
eta_deciles = deciles*-2.911;
gamma_deciles = deciles*.244;
for decile_idx = 1:length(deciles)
    for run = 1:nruns
        B = generative_model(seed, D, m, 'neighbors', {'powerlaw', 'powerlaw'}, ...
            [eta_deciles(decile_idx) gamma_deciles(decile_idx)]);
        b = zeros(n);
        b(B(:, 1)) = 1;
        b = b + b';
        synthetic_connectome_global_efficiency(decile_idx, run) = efficiency_bin(b);
        clear b B
        fprintf('Run %g of %g for decile %g\n', run, nruns, decile_idx);
    end
end
%% Specify neighbours generative network model function
function output = fcn_neighbors(A,K,D,m,eta,gam,epsilon)
% Note, this function only runs the power-law version of the neighbours
% generative model, as this is what we're collecting the output for.
step=1;
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
Fd = D.^eta;
Fk = K.^gam;
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);
% Initialise array to hold the networks at each step
Aall = zeros(m-mseed, n, n);
% Keep the first network
Aall(1,:,:) = A;
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
    Ff(uu,y) = Fd(uu,y).*(K(uu,y).^gam);
    Ff(y,uu) = Ff(uu,y)';
    Ff(vv,x) = Fd(vv,x).*(K(vv,x).^gam);
    Ff(x,vv) = Ff(vv,x)';
    Ff(A) = 0;
    P = Ff(indx);
    % Update the step and keep the networks!
    step=step+1;
    Aall(step,:,:) = A;
end
output = Aall;
end


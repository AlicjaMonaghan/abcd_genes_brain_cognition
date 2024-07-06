%% Generative variability for specific eta, gamma combinations
% Written by Danyal Akarca, University of Cambridge, 2023, and updated by
% Alicja Monaghan, University of Cambridge, 2024. The parcellation is
% defined as the Schaefer 100-node 17-network parcellation.
%% Set directories, load data, set paths
clear;clc;
cd('/Imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_variation_structural_generative_mechanisms_open/');
% Load the coordinates for the Schaefer 100-node parcellation, and
% calculate the Euclidean distance between points. 
schaefer100x17_1mm_info = load('data/schaefer100x17_1mm_info.mat');
schaefer100x17_1mm_info = schaefer100x17_1mm_info.schaefer100x17_1mm_info;
parcellation_coordinates = [schaefer100x17_1mm_info.x_mni,...
            schaefer100x17_1mm_info.y_mni, schaefer100x17_1mm_info.z_mni];
D = squareform(pdist(parcellation_coordinates));
% add bct
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
%% set key parameters
% number of connections
m = 331;
% Group 1: [-3.0610 0.2138]
% Group 2: [-2.9630 0.2144]
eta = -2.9630;
gamma = 0.2144;
% number of simulations
nrep = 1000;
%% compute current generative parameters
% form an even space based off the eta and gamma limits
x = eta; y = gamma;
params = [x y];
% set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'},{'powerlaw'}];
% take the nnode
nnode = size(D,1);
% number of connections
m = 331;
% load alicja's seed
load('/imaging/astle/users/da04/Postdoc/pgs_dissim/seed.mat');
Aseed = seed;
mseed = nnz(Aseed)/2;
% initialise
output = struct;
output.parameters = params;
output.networks = zeros(nrep,m);
output.probabiltiies = zeros(nrep,m-mseed);
%% run the generative model
fprintf('running network %g times %g connections...\n',nrep,m)
for rep = 1:nrep
    % display
    fprintf('running run %g of %g...\n',rep,nrep);
    % run the generative model
    % b = generative_model(Aseed,D,m,'matching',modelvar,params);
    epsilon = 1e-5;
    n = length(D);
    nparams = size(params,1);
    b = zeros(m,nparams);
    Pdev = zeros(m-mseed,nparams);
    Kseed = matching_ind(Aseed);
    Kseed = Kseed + Kseed';
    for iparam = 1:nparams
        eta = params(iparam,1);
        gam = params(iparam,2);
        outcome = fcn_matching(Aseed,Kseed,D,m,eta,gam,modelvar,epsilon);
        b(:,iparam) = outcome.connections;
        Pdev(:,iparam) = outcome.probabilities;
    end
    % store the output
    output.networks(rep,:) = b';
    output.probabilities(rep,:) = Pdev;
end
% save file - make sure to change network parameters according to which PGS
% group we're inspecting.
save('data/group_2_generative_model.mat','output','-v7.3');
%% matching
function output = fcn_matching(A,K,D,m,eta,gam,modelvar,epsilon)
        step = 1;
        K = K + epsilon;
        n = length(D);
        mseed = nnz(A)/2;
        mv1 = modelvar{1};
        mv2 = modelvar{2};
        switch mv1
            case 'powerlaw'
                Fd = D.^eta;
            case 'exponential'
                Fd = exp(eta*D);
        end
        switch mv2
            case 'powerlaw'
                Fk = K.^gam;
            case 'exponential'
                Fk = exp(gam*K);
        end
        Ff = Fd.*Fk.*~A;
        [u,v] = find(triu(ones(n),1));
        indx = (v - 1)*n + u;
        P = Ff(indx);
        Pdev = nan((n*(n-1)/2),m-mseed); % initialise
        for ii = (mseed + 1):m;

            C = [0; cumsum(P)];
            r = sum(rand*C(end) >= C);
            uu = u(r);
            vv = v(r);

            A(uu,vv) = 1;
            A(vv,uu) = 1;

            updateuu = find(A*A(:,uu));
            updateuu(updateuu == uu) = [];
            updateuu(updateuu == vv) = [];

            updatevv = find(A*A(:,vv));
            updatevv(updatevv == uu) = [];
            updatevv(updatevv == vv) = [];

            c1 = [A(:,uu)', A(uu,:)];
            for i = 1:length(updateuu)
                j = updateuu(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(uu) = 0;  use(uu+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(uu,j) = epsilon;
                    K(j,uu) = epsilon;
                else
                    K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,uu) = K(uu,j);
                end

            end

            c1 = [A(:,vv)', A(vv,:)];
            for i = 1:length(updatevv)
                j = updatevv(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(vv) = 0;  use(vv+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(vv,j) = epsilon;
                    K(j,vv) = epsilon;
                else
                    K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,vv) = K(vv,j);
                end
            end
            switch mv2
                case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            Ff = Fd.*Fk.*~A;
            P = Ff(indx); 
            Pdev(:,step) = P; % added
            step = step + 1;
        end
        b = find(triu(A,1));
        output = struct;
        output.connections = b;
        output.probabilities = mean(Pdev,1)';
 end
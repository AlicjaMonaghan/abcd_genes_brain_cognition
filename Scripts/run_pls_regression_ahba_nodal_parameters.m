function PLS_Output = run_pls_regression_ahba_nodal_parameters(subject)
% This script details a partial least squares (PLS) regression of 
% parameterised homophily-neighbours structural generative network model 
% (GNM)parameters for single ABCD individuals, one at a time, such as nodal
% wiring cost (PLS1) and nodal wiring value (PLS2), onto regional gene 
% expression values from the Allen Human Brain Atlas (AHBA). We restrict
% our analysis to the left hemisphere only, for which all 6 AHBA donors
% provided data for. We conduct permutation testing to draw null 
% distributions for the relationship between parameterised nodal parameters
% and AHBA gene arrays. The input for this function is a number between 1
% and 2153 which indexes the participant. Correspondence to Alicja
% Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk

% NOTE - The output for this function is the PLS_Output structure. In our
% study, we pulled the PLS outputs for each participant, as the PLSs would
% take too long to run in series. The
% linking_ahba_genetics_with_gnm_parameters_and_cognition.mat script takes
% the pulled PLS outputs for all participants, finds the genes which are
% highly predictive of parameterised nodal wiring costs or values across
% all participants. 

% NOTE - This code is based off of a version developed by Dr. Danyal
% Akarca.

% STEPS.
% 1. Setting up the work space.
% 2. Formatting data frames.
% 3. AHBA genes as predictors of parameterised nodal wiring costs - PLS1
% 4. AHBA genes as predictors of parameterised nodal wiring value - PLS2

%% PART 1 - Set up the workspace.
% Load up the directory for the parameterised nodal costs and values. <--
% SET THIS TO WHERE YOU SAVED THE DIRECTORY
cd('abcd_genomic_variation_structural_generative_mechanisms_open/');
% Loading the parameterised nodal wiring costs and values, respectively,
% for each participant.
parameterised_nodal_costs = load('data/neighbour_parameterised_nodal_costs.mat');
parameterised_nodal_costs = parameterised_nodal_costs.Neighbour_Parameterised_Nodal_Cost;
parameterised_nodal_value = load('data/neighbour_parameterised_nodal_value.mat');
parameterised_nodal_value = parameterised_nodal_value.Neighbour_Parameterised_Nodal_Value;
fprintf('We are processing participant number %d!\n', subject);
% Now load up the regional AHBA gene arrays. Genes were filtered according
% to those that had associated RNA-sequencing data from 2 AHBA donors,
% using the abagen toolbox (Markello et al., 2021; Arnatkeviciute et al.,
% 2019; Hawrylycz et al., 2012).
ahba_gene_expression = table2array(readtable('data/AHBA_expression_rnaseq_schaefer100_cleaned.csv', 'VariableNamingRule','preserve'));

%% PART 2 - Format Data Frames.
% We shall only focus on the left hemisphere, which all 6 AHBA participants
% have data for. Therefore, the dimensions of the AHBA_Gene_Expression_LH
% dataset will be 50 x 12432, whilst the dimensions of the
% Parameterised_Nodal_Cost_LH and Parameterised_Nodal_Value_LH datasets
% shall be nsub x 50 x 50, respectively. Clear up variables as we go
% along to save space. First, we need to select the parameterised nodal
% costs and value for the specific participant being processed, and only
% the first 50 nodes. 
ahba_gene_expression_lh = ahba_gene_expression(1:50,:);
parameterised_nodal_costs_lh = squeeze(parameterised_nodal_costs(subject,1:50,1:50));
parameterised_nodal_value_lh = squeeze(parameterised_nodal_value(subject,1:50,1:50));

% For the parameterised nodal costs and values in the left hemisphere, sum
% up the columns, excluding the diagonal, so that we find the total
% parameterised costs and values from each node to each other node. First,
% we need to initialise output arrays.
nnode = 50;
parameterised_nodal_costs_lh_summed = zeros(nnode,1);
parameterised_nodal_value_lh_summed = zeros(nnode,1);
% Remove the diagonal
parameterised_nodal_costs_lh = squeeze(parameterised_nodal_costs_lh) - diag(diag(squeeze(parameterised_nodal_costs_lh)));
% And sum the columns for each participant
parameterised_nodal_costs_lh_summed(:,1) = squeeze(sum(parameterised_nodal_costs_lh));
% Repeat the same process for parameterised nodal value.
parameterised_nodal_value_lh = squeeze(parameterised_nodal_value_lh) - diag(diag(squeeze(parameterised_nodal_value_lh)));
parameterised_nodal_value_lh_summed(:,1) = squeeze(sum(parameterised_nodal_value_lh));

% Clear the old variables to save space
clear ahba_gene_expression parameterised_nodal_value parameterised_nodal_costs parameterised_nodal_value_lh parameterised_nodal_costs_lh AHBA_Missing_Row

%% PART 3 - Conduct Nodal Costs Partial Least Squares Regression!
% For each participant, we shall conduct two PLSs i.e. PLS1 is
% Parameterised_Nodal_Cost_LH ~ AHBA_Gene_Expression_LH and PLS2 is
% Parameterised_Nodal_Value_LH ~ AHBA_Gene_Expression_LH. We shall conduct
% 10000 permutations to compare the test statistic with a null distribution.
% We shall also correct for multiple comparisons. First, initialise output
% arrays for the PLS, after setting the number of genes, PLS components, 
% and the permutations per gene. 
ngenes = width(ahba_gene_expression_lh);
ncomp = 3;
nperm = 10000;
Xloading = zeros(ngenes,ncomp);
Xscore = zeros(nnode,ncomp);
var = zeros(2,ncomp);
% Initialising outputs for the permutations.
Xloading_perm = zeros(nperm,ngenes,ncomp);
% Initialising corrected p values for gene loadings.
pcorr = zeros(ngenes,ncomp);

% Conduct the PLS of parameterised nodal costs onto AHBA gene expression. 
% Start the timer 
tic;
% Update user
fprintf('Starting to run the nodal cost PLS for participant %d\n', subject);
% Run the PLS for this participant. 
[Xloading,~,Xscore,~,~,var] = plsregress(ahba_gene_expression_lh,parameterised_nodal_costs_lh_summed,ncomp);
% Now set up the permutation for this participant.
for permutation = 1:nperm
    fprintf('Running permutation %g of %g for participant %d\n', permutation, nperm, subject);
    % Change the order of the response variable 
    Y = parameterised_nodal_costs_lh_summed(randperm(length(parameterised_nodal_costs_lh_summed)));
    % And run the actual permutation...
    [Xl,~,~,~,~,~] = plsregress(ahba_gene_expression_lh,Y,ncomp);
    % Keep the outputs from the permutation
    Xloading_perm(permutation,:,:) = Xl;
    % Clear variables to save space
    clear Y Xl 
end
% Calculate the p vector for this participant for all components
for comp = 1:ncomp
    % Select the null matrix for this component
    null = squeeze(squeeze(Xloading_perm(:,:,comp)));
    % Now loop across each of the 15792 genes
    for gene = 1:ngenes
        % Gather the observed X loading for this gene
        r = Xloading(gene,comp);
        % Gather the associated null value
        z = null(:,gene)';
        % Find the p vector!
        pcorr(gene,comp) = 1-(sum(z<r)/nperm);
    end
    % Clear the null distribution and other variables to save space
    clear null r z
end
Participant_Processing_Time = toc;
% Update the user
fprintf('Completed observed and permutated PLS for nodal cost for participant %d in %g seconds!\n', ...
    subject, round(Participant_Processing_Time,2));

% Assign the outputs to a structure for nodal value.
PLS_Output = struct();
PLS_Output.Nodal_Cost_PLS = struct();
PLS_Output.Nodal_Cost_PLS.pcorr = pcorr;
PLS_Output.Nodal_Cost_PLS.var = var;
PLS_Output.Nodal_Cost_PLS.Xloading = Xloading;
PLS_Output.Nodal_Cost_PLS.Xscore = Xscore;
PLS_Output.Nodal_Cost_PLS.Permutation_Xloading = Xloading_perm; 
% Clear the above variables to save space
clear Xloading Xscore var Xloading_perm pcorr PLS_Output.Nodal_Costs_PLS
%% PART 4 - Conduct Nodal Value Partial Least Squares Regression!
% This is the second PLS we shall run for each participant. First, 
% initialise output arrays, where the number of genes, permutations etc 
% remain the same from Part 3. 
Xloading = zeros(ngenes,ncomp);
Xscore = zeros(nnode,ncomp);
var = zeros(2,ncomp);
% Initialising outputs for the permutations.
Xloading_perm = zeros(nperm,ngenes,ncomp);
% Initialising corrected p values for gene loadings.
pcorr = zeros(ngenes,ncomp);

% Conduct the PLS of parameterised nodal value onto AHBA gene expression. 
% Start the timer 
tic;
% Update user
fprintf('Starting to run the nodal value PLS for participant %d\n',subject);
% Run the PLS for this participant. 
[Xloading,~,Xscore,~,~,var] = plsregress(ahba_gene_expression_lh,parameterised_nodal_value_lh_summed,ncomp);
% Now set up the permutation for this participant.
for permutation = 1:nperm
    fprintf('Running permutation %g of %g for participant %d\n',permutation,nperm,subject);
    % Change the order of the response variable 
    Y = parameterised_nodal_value_lh_summed(randperm(length(parameterised_nodal_value_lh_summed)));
    % And run the actual permutation...
    [Xl,~,~,~,~,~] = plsregress(ahba_gene_expression_lh,Y,ncomp);
    % Keep the outputs from the permutation
    Xloading_perm(permutation,:,:) = Xl;
    % Clear variables to save space
    clear Y Xl 
end
% Calculate the p vector for this participant for all components
for comp = 1:ncomp
    % Select the null matrix for this component
    null = squeeze(squeeze(Xloading_perm(:,:,comp)));
    % Now loop across each of the 15792 genes
    for gene = 1:ngenes
        % Gather the observed X loading for this gene
        r = Xloading(gene,comp);
        % Gather the associated null value
        z = null(:,gene)';
        % Find the p vector!
        pcorr(gene,comp) = 1-(sum(z<r)/nperm);
    end
    % Clear the null distribution and other variables to save space
    clear null r z
end
Participant_Processing_Time = toc;
% Update the user
fprintf('Completed observed and permutated PLS for nodal value for participant %d in %g seconds!\n', ...
    subject, round(Participant_Processing_Time,2));

% Assign the outputs to a structure for nodal costs.
PLS_Output.Nodal_Value_PLS = struct();
PLS_Output.Nodal_Value_PLS.pcorr = pcorr;
PLS_Output.Nodal_Value_PLS.var = var;
PLS_Output.Nodal_Value_PLS.Xloading = Xloading;
PLS_Output.Nodal_Value_PLS.Xscore = Xscore;
PLS_Output.Nodal_Value_PLS.Permutation_Xloading = Xloading_perm; 

% Clear the above variables to save space
clear Xloading Xscore var Xloading_perm pcorr PLS_Output
% And update the user
fprintf('The two PLSs for participant %d are complete!\n', subject);
% And save!
save(sprintf('data/ahba_gnm_nodal_pls_participant_%d.mat',subject),"PLS_Output");
end

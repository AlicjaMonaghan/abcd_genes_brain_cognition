% This script details linking nodal wiring parameters (eta and gamma), from
% generative network models (GNM), with cognitive ability scores and 
% polygenic scores (PGS) for cognitive ability, in the Adolescent Brain 
% Cognitive Development (ABCD) data set. We include sex, in-scanner motion,
% scanner site and age as covariates. First, we assess the relationship 
% between GNM parameters and cognitive ability using a partial least 
% squares regression (PLS), with the former as predictors and covariates, 
% and the latter as the response. Results for all 3 PLSs are reported in
% Supplementary Table 8.  We then sort
% the genes by decreasing X loading in order to submit them for gene
% ontology analysis, using the abcd_gene_ontology_analysis.R script. All
% correspondence to Alicja Monaghan, alicja.monaghan@mrc-cbu.cam.ac.uk

% STEPS:
% 1. Preparing the work space.
% 2. Exploratory univariate analyses between different levels of
% explanation i.e. genetic (PGSs for cognitive ability), brain (GNM
% parameters eta and gamma), and behaviour (general g factor of cognition).
% 3. PLS with GNM parameters (eta and gamma) as predictors of cognitive
% ability. 
% 4. PLS with GNM parameters as predictors of PGSs for cognitive ability.
% 5. PLS with GNM parameters and PGSs for cognitive ability as predictors
% of cognitive ability.

%% PART 1 - Set up the workspace.
% Load up the ABCD individual eta and gamma parameters with polygenic
% scores, and format! Convert the variables to numeric as appropriate.
clear;clc;
cd('//cbsu/data/Imaging/projects/external/abcd/analyses/Alicja/abcd_genomic_variation_structural_generative_mechanisms_open/');
% Load GNM parameters (eta and gamma), polygenic scores for general
% cognitive ability, cognitive ability scores, and covariates.
parameters_and_covariates = readtable('data/gnm_parameters_pgs_covariates_anonymised.csv');

% Extract the individual variables from the data.
eta = parameters_and_covariates.eta;
gamma = parameters_and_covariates.gamma;
pgs = parameters_and_covariates.pgs;
cognitive_ability = parameters_and_covariates.pc1;
age = parameters_and_covariates.age;
scanner_site = char(parameters_and_covariates.siteid);
mean_fwd = parameters_and_covariates.meanFWD;
sex = parameters_and_covariates.sex;

% Find how many participants have genetic, neuroimaging, and cognitive
% data!
nsub_genetic_neuro_cog = length(eta);

% Load the AHBA gene expression arrays.
ahba_gene_expression = readtable('data/AHBA_expression_rnaseq_schaefer100_cleaned.csv',"VariableNamingRule","preserve");
% Retrieve the probe names associated with each gene
ahba_gene_names = ahba_gene_expression.Properties.VariableNames;
% Remove the label column.
ahba_gene_names(:,1) = [];

%% PART 2 - Exploratory Univariate Analyses Between Different Levels of Explanation %%
% Correlate eta and gamma.
[rho,pval] = corr(eta,gamma,'Type','spearman');
fprintf('Spearman-rank correlation between eta and gamma is r(1459) = %.03d, p = %.03d\n',rho,pval);
% Correlate eta with cognitive ability.
[rho,pval] = corr(eta,cognitive_ability,'Type','spearman');
fprintf('Spearman-rank correlation between eta and cognitive ability is r(1459) = %.03d, p = %.03d\n',rho,pval);
% Correlate gamma with cognitive ability.
[rho,pval] = corr(gamma,cognitive_ability,'Type','spearman');
fprintf('Spearman-rank correlation between gamma and cognitive ability is r(1459) = %.03d, p = %.03d\n',rho,pval);
% Correlate eta with PGSs for cognitive ability.
[rho, pval] = corr(eta, pgs, 'Type', 'spearman');
fprintf('Spearman-rank correlation between eta and PGSs for cognitive ability is r(1459) = %.03d, p = %.03d\n',rho,pval);
% Correlate gamma with PGSs for cognitive ability.
[rho, pval] = corr(gamma, pgs, 'Type', 'spearman');
fprintf('Spearman-rank correlation between gamma and PGSs for cognitive ability is r(1459) = %.03d, p = %.03d\n',rho,pval);

%% PART 3 - Examine Relationship Between GNM Parameters and Cognitive Ability Using PLS %%
% To examine the relationship between PGSs for cognitive ability and GNM
% parameters, we shall conduct a PLS of GNM wiring parameters as predictors
% for cognitive ability. We shall include 4 covariates: age, sex, scanner 
% site and in-scanner motion. To draw a null distribution for significance 
% testing, we shall conduct 10,000 permutations. Finally, to derive mean 
% loadings and confidence intervals, we shall bootstrap the permuted X 
% loadings 10,000 times. As an additional robustness check, we shall 
% compare the root mean square error (RMSE) of the observed model compared 
% to the 10,000 null (permuted) models. 

% First, normalise all continuous predictor and response variables, then
% set these as X and Y, respectively.
X = [normalize([eta gamma])];
% Normalise the outcome variable.
Y = normalize(cognitive_ability);
% Normalise the continuous covariates.
Covariates = [normalize([age mean_fwd]) scanner_site sex];
% Now run the PLS, with 10-fold cross-validation, and the default number of
% components. According to the plsregress vigenette, the default ncomp is
% the smaller of size(X,1) – 1 and size(X,2), so in our case is 2!
ncomp = size(X,2);
[XL_Observed,YL_Observed,XS_Observed,YS_Observed,Beta,PCTVAR_Observed,MSE_Observed_Eta_Gamma_Predicting_Cognition,Stats_Observed] = plsregress(X,Y,ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,'KFold',10));
% Get the root mean square error (RMSE).
RMSE_Observed_Eta_Gamma_Predicting_Cognition = sqrt(MSE_Observed_Eta_Gamma_Predicting_Cognition(2,end));

% And find the correlation between the component scores for the observed
% PLS for each component after finding how many components were extracted.
% Adjust for the covariates too through partial correlation.
ncomp = size(XL_Observed,2);
Observed_Corr = zeros(ncomp,1);
for comp = 1:ncomp
    Observed_Corr(comp) = partialcorr(XS_Observed(:,comp), YS_Observed(:,comp), Covariates);
end
% Now initiate the permuted output variables for the number of PLS
% components extracted (2). nperm will be the same for all PLSs, and
% therefore is specified here once. 
nperm = 10000;
Permuted_XS = zeros(size(XS_Observed,1),size(XS_Observed,2),nperm);
Permuted_YS = zeros(size(YS_Observed,1),size(YS_Observed,2),nperm);
Permuted_Xloadings = zeros(size(XL_Observed,1),size(XL_Observed,2),nperm);
Permuted_RMSE_Eta_Gamma_Predicting_Cognition = zeros(nperm,1); 
% Run the permutations, after setting the seed!
rng('default');
for permutation = 1:nperm
    % Shuffle and randomise the predictors
    Randomised_Predictors = X(randperm(length(X)),:);
    % Run the permuted PLS, and extract the required parameters! Note that 
    % we normalised the outcome variable at the start. 
    [XL_Perm,~,XS_Perm,YS_Perm,~,~,MSE_Perm] = plsregress(Randomised_Predictors,Y,ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,'KFold',10));
    % And append the XS, and YS to the permuted output array, as well as 
    % the permuted X loadings (for bootstrapping). 
    Permuted_XS(:,:,permutation) = XS_Perm;
    Permuted_YS(:,:,permutation) = YS_Perm;
    Permuted_Xloadings(:,:,permutation) = XL_Perm;
    % Get the RMSE and append to the output array
    Permuted_RMSE_Eta_Gamma_Predicting_Cognition(permutation,1) = sqrt(MSE_Perm(2,end));
    clear XS_Perm YS_Perm XL_Perm MSE_Perm
end

% Now correlate the permuted XS and permuted YS vectors across components, 
% after initialising the Permuted_Corr variable, and adjusting for 
% covariates. We will then compare the observed XS-YS correlations with 
% their permuted counterparts for each component, after initialising an 
% output variable. 
Permuted_Corr = zeros(ncomp,nperm);
pcorr = zeros(ncomp,1);
for comp = 1:ncomp
    % Correlate the permuted X and Y scores.
    Permuted_Correlation_Matrix = partialcorr(squeeze(Permuted_XS(:,comp,:)),squeeze(Permuted_YS(:,comp,:)), Covariates);
    % Select the first column of this matrix and assign to the output.
    Permuted_Corr(comp,:) = Permuted_Correlation_Matrix(:,1);
    % Extract the null distribution for this component.
    null_distribution = Permuted_Corr(comp,:);
    % Extract the observed PLS XS-YS correlation.
    Observed = Observed_Corr(comp,:);
    % And compare the extracted and observed XS-YS correlations.
    pcorr(comp,:) = 1 - (sum(null_distribution<Observed)/nperm);
end
% Now find test the RMSE distribition!
pcorr_RMSE_Eta_Gamma_Predicting_Cognition = 1 - (sum(Permuted_RMSE_Eta_Gamma_Predicting_Cognition<RMSE_Observed_Eta_Gamma_Predicting_Cognition)/nperm);

% For each component, we shall bootstrap the PLS loadings to extract the 
% mean and confidence intervals of each predictor. Initialise the output
% variables first!
npred = size(X,2);
% For each Cognition_Onto_Eta_Gamma_PLS_CI array, the first of the third dimensions
% represents the lower limit, whilst the second of the third dimension
% represents the upper limit. 
Cognition_Onto_Eta_Gamma_PLS_CI_Bootstrapped_Xloadings = zeros(npred,ncomp,2);
Cognition_Onto_Eta_Gamma_PLS_CI_Bootstrapped_Yloadings = zeros(ncomp,2);
% We shall now conduct bootstrapping of the original PLS loadings to
% extract mean loadings and confidence intervals. Draw n samples with
% replacement and repeat the PLS. First, set the number of bootstraps we
% shall conduct, and initialise the output variable for the X and Y 
% loadings. nboot will be the same for all PLSs, and therefore is specified
% here once. 
nboot = 10000;
XL_Bootstrapped_Array = zeros(size(XL_Observed,1),size(XL_Observed,2),nboot);
YL_Bootstrapped_Array = zeros(size(YL_Observed,1),size(YL_Observed,2),nboot);
Corrected_XL_Bootstrapped = zeros(size(XL_Observed,1),size(XL_Observed,2));
Corrected_YL_Bootstrapped = zeros(size(YL_Observed,1),size(YL_Observed,2));
% Now conduct the bootstrapped PLSs!
for boot = 1:nboot
    % Randomly sample nsub participants with replacement and normalize.
    Randomly_Sampled_Data = normalize(datasample([eta gamma cognitive_ability], nsub_genetic_neuro_cog));
    % Repeat the PLS and keep the X loadings.
    [XL_Bootstrapped, YL_Bootstrapped] = plsregress(Randomly_Sampled_Data(:,1:npred),Randomly_Sampled_Data(:,end),ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,"KFold",10));
    % Assign to the output variables. Make sure to account for sign 
    % flipping using the Procrustes correction before assigning to the
    % output for each component!. The target for the correction will be the
    % observed X and Y loadings. 
    for comp = 1:ncomp
            Corrected_XL_Bootstrapped(:,comp) = rotatefactors(XL_Bootstrapped(:,comp),'Method','procrustes','Target',XL_Observed(:,comp));
            XL_Bootstrapped_Array(:,comp,boot) = squeeze(Corrected_XL_Bootstrapped(:,comp));
            Corrected_YL_Bootstrapped(:,comp) = rotatefactors(YL_Bootstrapped(:,comp),'Method','procrustes','Target',YL_Observed(:,comp));
            YL_Bootstrapped_Array(:,comp,boot) = squeeze(Corrected_YL_Bootstrapped(:,comp));
    end
    % Clear the XL_Bootstrapped variable and repeat
    clear Corrected_XL_Bootstrapped YL_Bootstrapped
end
YL_Bootstrapped_Array = squeeze(YL_Bootstrapped_Array);
% Now calculate the mean X and Y loadings for each predictor and component.
Cognition_Onto_Eta_Gamma_PLS_Mean_X = mean(XL_Bootstrapped_Array,3);
Cognition_Onto_Eta_Gamma_PLS_Mean_Y = mean(YL_Bootstrapped_Array,2);
% Find the confidence interval for the mean bootstrapped X and Y loadings 
% for each predictor and component.
for comp = 1:ncomp
    for pred = 1:npred
        % Extract the X and Y loadings for this predictor and component.
        Bootstrapped_XL_Pred_Comp = squeeze(XL_Bootstrapped_Array(pred,comp,:));
        Bootstrapped_YL_Comp = squeeze(YL_Bootstrapped_Array(comp,:));
        % Calculate the standard error for the X and Y loadings!
        SEM_X = std(Bootstrapped_XL_Pred_Comp)/sqrt(length(Bootstrapped_XL_Pred_Comp));
        SEM_Y = std(Bootstrapped_YL_Comp)/sqrt(length(Bootstrapped_YL_Comp));
        % Calculate the t-score for a 95% confidence interval.
        ts_X = tinv([0.025 0.975], length(Bootstrapped_XL_Pred_Comp)-1);
        ts_Y = tinv([0.025 0.975], length(Bootstrapped_YL_Comp)-1);
        % And calculate the confidence intervals for each!
        CI_X = mean(Bootstrapped_XL_Pred_Comp) + ts_X*SEM_X;
        CI_Y = mean(Bootstrapped_YL_Comp) + ts_Y*SEM_Y;
        % Assign to the output variable
        Cognition_Onto_Eta_Gamma_PLS_CI_Bootstrapped_Xloadings(pred,comp,:) = CI_X;
        Cognition_Onto_Eta_Gamma_PLS_CI_Bootstrapped_Yloadings(comp,:) = CI_Y;
    end 
end

% Clear up intermediate variables!
clear Bootstrapped_XL_Pred_Comp Bootstrapped_YL_Comp SEM_X SEM_Y ts_X ts_Y 
clear CI_X CI_Y 

% We shall now find which of the parameters, eta or gamma, is particularly
% important in the PLS, by calculating the variable importance in 
% projection (VIP) scores. Note that the following uses code from the 
% plsregress MATLAB vignette.

% Calculate the normalised beta weights from the observed PLS.
W0 = Stats_Observed.W ./ sqrt(sum(Stats_Observed.W.^2,1));
% Calculate the VIP scores for the two components.
p = size(XL_Observed,1);
sumSq = sum(XS_Observed.^2,1).*sum(YL_Observed.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
% And find variables with a VIP score greater than 1. This gives us an
% index of 2 for the predictor variables, indicating gamma.
indVIP = find(vipScore >= 1);

%% PART 4 - Linking GNM Parameters and PGS for Cognitive Ability
% Now, to investigate whether eta and gamma together are reflected in the
% polygenic score, we shall conduct a PLS of eta and gamma as the two
% predictors, and PGS as the response variable. Note that the X vector and
% covariates have already been normalised from the section above, so we
% only need to normalise the PGS here.
Y = normalize(pgs);
% Run the (observed) PLS, with the default number of components and 10-fold 
% cross-validation.
ncomp = size(X,2);
[XL_Observed,YL_Observed,XS_Observed,YS_Observed,~,PCTVAR_Observed,MSE_Observed_Eta_Gamma_Predicting_PGS,Stats_Observed] = plsregress(X,Y,ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,'KFold',10));
% Find the RMSE
RMSE_Eta_Gamma_Predicting_PGS = sqrt(MSE_Observed_Eta_Gamma_Predicting_PGS(2,end));

% And find the correlation between the component scores for the observed
% PLS for each component, after finding how many components were extracted.
ncomp = size(XL_Observed,2);
Observed_Corr = zeros(ncomp,1);
for comp = 1:ncomp
    Observed_Corr(comp) = partialcorr(XS_Observed(:,comp), YS_Observed(:,comp), Covariates);
end

% Now initiate the permuted output variables for each component.
Permuted_XS = zeros(size(XS_Observed,1),size(XS_Observed,2),nperm);
Permuted_YS = zeros(size(YS_Observed,1),size(YS_Observed,2),nperm);
Permuted_RMSE_Eta_Gamma_Predicting_PGS = zeros(nperm,1); 

% Run the permutations, after setting the seed!
rng('default');
for permutation = 1:nperm
    % Shuffle and randomise the predictors
    Randomised_Predictors = X(randperm(length(X)),:);
    % Run the permuted PLS, and extract the required parameters! Note that 
    % we normalised the outcome variable at the start. 
    [~,~,XS_Perm,YS_Perm,~,~,MSE_Perm] = plsregress(Randomised_Predictors,Y,ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,"KFold",10));
    % And append the XS and YS variables to the permuted output array.
    Permuted_XS(:,:,permutation) = XS_Perm;
    Permuted_YS(:,:,permutation) = YS_Perm;
    % And append the RMSE
    Permuted_RMSE_Eta_Gamma_Predicting_PGS(permutation,1) = sqrt(MSE_Perm(2,end));
    % And clear the intermediate variables!
    clear XS_Perm YS_Perm MSE_Perm
end
% Now find test the RMSE distribition!
pcorr_RMSE_Eta_Gamma_Predicting_PGS = 1 - (sum(Permuted_RMSE_Eta_Gamma_Predicting_PGS<RMSE_Eta_Gamma_Predicting_PGS)/nperm);

% Now correlate the permuted XS and permuted YS vectors across components, 
% after initialising the Permuted_PCorr variable. We will then compare the
% observed XS-YS correlations with their permuted counterparts for each
% component, after initialising an output variable. 
Permuted_Corr = zeros(ncomp,nperm);
pcorr = zeros(ncomp,1);
for comp = 1:ncomp
    % Correlate the permuted X and Y scores.
    Permuted_Correlation_Matrix = partialcorr(squeeze(Permuted_XS(:,comp,:)),squeeze(Permuted_YS(:,comp,:)), Covariates);
    % Select the first column of this matrix and assign to the output.
    Permuted_Corr(comp,:) = Permuted_Correlation_Matrix(:,1);
    % Extract the null distribution for this component.
    null_distribution = Permuted_Corr(comp,:);
    % Extract the observed PLS XS-YS correlation.
    Observed = Observed_Corr(comp,:);
    % And compare the extracted and observed XS-YS correlations.
    pcorr(comp,:) = 1 - (sum(null_distribution<Observed)/nperm);
end

% For each component, bootstrap the PLS loadings to extract the mean and 
% confidence intervals of each predictor. npred is the same as Part 2. 
Eta_Gamma_PGS_PLS_CI_X = zeros(npred,ncomp,2);
Eta_Gamma_PGS_PLS_CI_Y = zeros(ncomp,2);

% Bootstrap the original PLS loadings to extract the mean and confidence 
% intervals. Draw n samples with replacement and repeat the PLS. 
XL_Bootstrapped_Array = zeros(size(XL_Observed,1),size(XL_Observed,2),nboot);
YL_Bootstrapped_Array = zeros(size(YL_Observed,1),size(YL_Observed,2),nboot);

% Now conduct the bootstrapped PLSs!
for boot = 1:nboot
    % Randomly sample nsub participants with replacement and normalize.
    Randomly_Sampled_Data = normalize(datasample([X Y], nsub_genetic_neuro_cog));
    % Repeat the PLS and keep the X loadings.
    [XL_Bootstrapped, YL_Bootstrapped] = plsregress(Randomly_Sampled_Data(:,1:npred),Randomly_Sampled_Data(:,end),ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,"KFold",10));
    % Correct for sign flipping using the Procrustes correction applied
    % across each component. 
    for comp = 1:ncomp
        Corrected_XL_Bootstrapped(:,comp) = rotatefactors(XL_Bootstrapped(:,comp),'Method','procrustes','Target',XL_Observed(:,comp));
        XL_Bootstrapped_Array(:,comp,boot) = squeeze(Corrected_XL_Bootstrapped(:,comp));
        Corrected_YL_Bootstrapped(:,comp) = rotatefactors(YL_Bootstrapped(:,comp),'Method','procrustes','Target',YL_Observed(:,comp));
        YL_Bootstrapped_Array(:,comp,boot) = squeeze(Corrected_YL_Bootstrapped(:,comp));
    end
    % Clear the intermediate variables and repeat
    clear Corrected_XL_Bootstrapped Corrected_YL_Bootstrapped
end

YL_Bootstrapped_Array = squeeze(YL_Bootstrapped_Array);
% Now calculate the mean X and Y loadings for each predictor and component.
Eta_Gamma_PGS_PLS_Mean_X = mean(XL_Bootstrapped_Array,3);
Eta_Gamma_PGS_PLS_Mean_Y = mean(YL_Bootstrapped_Array,2);
% Find the confidence interval for the mean bootstrapped X and Y loadings 
% for each predictor and component.
for comp = 1:ncomp
    for pred = 1:npred
        % Extract the X and Y loadings for this predictor and component.
        Bootstrapped_XL_Pred_Comp = squeeze(XL_Bootstrapped_Array(pred,comp,:));
        Bootstrapped_YL_Comp = squeeze(YL_Bootstrapped_Array(comp,:));
        % Calculate the standard error for each type of loading
        SEM_X = std(Bootstrapped_XL_Pred_Comp)/sqrt(length(Bootstrapped_XL_Pred_Comp));
        SEM_Y = std(Bootstrapped_YL_Comp)/sqrt(length(Bootstrapped_YL_Comp));
        % Calculate the t-score for a 95% confidence interval.
        ts_X = tinv([0.025 0.975], length(Bootstrapped_XL_Pred_Comp)-1);
        ts_Y = tinv([0.025 0.975], length(Bootstrapped_YL_Comp)-1);
        % And calculate the confidence interval!
        CI_X = mean(Bootstrapped_XL_Pred_Comp) + ts_X*SEM_X;
        CI_Y = mean(Bootstrapped_YL_Comp) + ts_Y*SEM_Y;
        % Assign to the output variable
        Eta_Gamma_PGS_PLS_CI_X(pred,comp,:) = CI_X;
        Eta_Gamma_PGS_PLS_CI_Y(comp,:) = CI_Y;
    end
end
 
% Clear up intermediate variables!
clear Bootstrapped_XL_Pred_Comp Bootstrapped_YL_Comp SEM_X SEM_Y ts_X ts_Y 
clear CI_X CI_Y 
 
% Calculate the normalised beta weights from the observed PLS.
W0 = Stats_Observed.W ./ sqrt(sum(Stats_Observed.W.^2,1));
% Calculate the VIP scores for the two components.
p = size(XL_Observed,1);
sumSq = sum(XS_Observed.^2,1).*sum(YL_Observed.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
% And find variables with a VIP score greater than 1. This gives us an
% index of 1 for the predictor variables, indicating eta. 
indVIP = find(vipScore >= 1);

%% PART 5 - Examine Relationship Between PGS, GNM Parameters, and Cognitive Ability!
% We shall now conduct a PLS of ABCD PGSs for cognitive ability and GNM 
% parameters as predictors, and cognitive ability as the response,
% correcting for covariates. 

% First, initialise these variables, after normalizing!
X = [normalize([eta gamma pgs])];
Y = normalize(cognitive_ability);
% Run the (observed) PLS with 10-fold cross-validation.
ncomp = size(X,2);
[XL_Observed,YL_Observed,XS_Observed,YS_Observed,~,PCTVAR_Observed,MSE_Observed_Eta_Gamma_PGS_Predicting_Cognition,Stats_Observed] = plsregress(X,Y,ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,"KFold",10));
% Extract the RMSE
RMSE_Eta_Gamma_PGS_Predicting_Cognition = sqrt(MSE_Observed_Eta_Gamma_PGS_Predicting_Cognition(2,end));
% And find the correlation between the component scores for the observed
% PLS for each component. Correct for covariates, which have already been
% normalised. 
ncomp = size(XL_Observed,2);
Observed_Corr = zeros(ncomp,1);
for comp = 1:ncomp
    Observed_Corr(comp) = partialcorr(XS_Observed(:,comp), YS_Observed(:,comp), Covariates);
end
% Now initiate the permuted output variables. Keep the permuted X loadings 
% too for later bootstrapping analyses. 
Permuted_XS = zeros(size(XS_Observed,1),size(XS_Observed,2),nperm);
Permuted_YS = zeros(size(YS_Observed,1),size(YS_Observed,2),nperm);
Permuted_Xloadings = zeros(size(XL_Observed,1),size(XL_Observed,2),nperm);
Perm_RMSE_Eta_Gamma_PGS_Predicting_Cognition = zeros(nperm,1);
% Run the permutations, after setting the seed.
rng('default');
for permutation = 1:nperm
    % Shuffle the predictors
    Randomised_Predictors = X(randperm(length(X)),:);
    % Run the permuted PLS and extract the XS and YS scores! Note that the 
    % outcome variable Y was already normalized in the previous step. 
    [XL_perm,~,XS_perm,YS_perm,~,~,MSE_Perm_Eta_Gamma_PGS_Predicting_Cognition] = plsregress(Randomised_Predictors,Y,ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,"KFold",10));
    % And append the XS and YS variables to the permuted output array, as
    % well as the permuted X loadings and RMSE.
    Permuted_XS(:,:,permutation) = XS_perm;
    Permuted_YS(:,:,permutation) = YS_perm;
    Permuted_Xloadings(:,:,permutation) = XL_perm;
    Perm_RMSE_Eta_Gamma_PGS_Predicting_Cognition(permutation) = sqrt(MSE_Perm_Eta_Gamma_PGS_Predicting_Cognition(2,end));
    % And clear intermediate variables
    clear XS_perm YS_perm XL_perm MSE_Perm_Eta_Gamma_PGS_Predicting_Cognition
end
% Test for a significant difference between the observed and permuted RMSE.
pcorr_RMSE_Eta_Gamma_PGS_Predicting_Cognition = 1 - (sum(Perm_RMSE_Eta_Gamma_PGS_Predicting_Cognition<RMSE_Eta_Gamma_PGS_Predicting_Cognition)/nperm);
% Now correlate the permuted XS and permuted YS vectors across components,
% after initialising the Permuted_PCorr variable, and correct for
% covariates through partial correlation. 
Permuted_Corr = zeros(ncomp,nperm);
pcorr = zeros(ncomp,1);
for comp = 1:ncomp
    % Correlate the permuted X and Y scores.
    Permuted_Correlation_Matrix = partialcorr(squeeze(Permuted_XS(:,comp,:)),squeeze(Permuted_YS(:,comp,:)), Covariates);
    % Select the first column of this matrix and assign to the output.
    Permuted_Corr(comp,:) = Permuted_Correlation_Matrix(:,1);
    % Extract the null distribution for this component.
    null_distribution = Permuted_Corr(comp,:);
    % Extract the observed PLS XS-YS correlation.
    Observed = Observed_Corr(comp,:);
    % And compare the extracted and observed XS-YS correlations.
    pcorr(comp,:) = 1 - (sum(null_distribution<Observed)/nperm);
end
% For each component, bootstrap the PLS loadings to extract the mean and
% confidence intervals of each predictor. Initialise the output variables.
npred = size(X,2);
Cognition_Onto_Eta_Gamma_PGS_PLS_CI_X = zeros(npred,ncomp,2);
Cognition_Onto_Eta_Gamma_PGS_PLS_CI_Y = zeros(ncomp,2);
% We shall now conduct bootstrapping of the original PLS loadings to
% extract mean loadings and confidence intervals.
XL_Bootstrapped_Array = zeros(size(XL_Observed,1),size(XL_Observed,2),nboot);
YL_Bootstrapped_Array = zeros(size(YL_Observed,1),size(YL_Observed,2),nboot);
Corrected_XL_Bootstrapped = zeros(size(XL_Observed,1),1);
Corrected_YL_Bootstrapped = zeros(size(1));
% Now conduct the bootstrapped PLSs!
for boot = 1:nboot
    % Randomly sample nsub participants with replacement and normalize.
    Randomly_Sampled_Data = normalize(datasample([X Y], nsub_genetic_neuro_cog));
    % Repeat the PLS and keep the X loadings.
    [XL_Bootstrapped, YL_Bootstrapped] = plsregress(Randomly_Sampled_Data(:,1:npred),Randomly_Sampled_Data(:,end),ncomp,"cv",cvpartition(nsub_genetic_neuro_cog,"KFold",10));
    % Assign to the output variables. Make sure to account for sign 
    % flipping using the Procrustes correction before assigning to the
    % output for each component!. The target for the correction will be the
    % observed X and Y loadings. 
    for comp = 1:ncomp
            Corrected_XL_Bootstrapped(:,comp) = rotatefactors(XL_Bootstrapped(:,comp),'Method','procrustes','Target',XL_Observed(:,comp));
            XL_Bootstrapped_Array(:,comp,boot) = squeeze(Corrected_XL_Bootstrapped(:,comp));
            Corrected_YL_Bootstrapped(:,comp) = rotatefactors(YL_Bootstrapped(:,comp),'Method','procrustes','Target',YL_Observed(:,comp));
            YL_Bootstrapped_Array(:,comp,boot) = squeeze(Corrected_YL_Bootstrapped(:,comp));
    end
    % Clear the intermediate variables and repeat
    clear Corrected_XL_Bootstrapped Corrected_YL_Bootstrapped
end
YL_Bootstrapped_Array = squeeze(YL_Bootstrapped_Array);
% Now calculate the mean X loadings for each predictor and component.
Cognition_Onto_Eta_Gamma_PGS_PLS_Mean_X = mean(XL_Bootstrapped_Array,3);
Cognition_Onto_Eta_Gamma_PGS_PLS_Mean_Y = mean(YL_Bootstrapped_Array,2);
% Find the confidence interval for the mean bootstrapped X loadings for
% each predictor and component.
for comp = 1:ncomp
    for pred = 1:npred
        % Extract the X and Y loadings for this predictor and component.
        Bootstrapped_XL_Pred_Comp = squeeze(XL_Bootstrapped_Array(pred,comp,:));
        Bootstrapped_YL_Comp = squeeze(YL_Bootstrapped_Array(comp,:));
        % Calculate the standard error
        SEM_X = std(Bootstrapped_XL_Pred_Comp)/sqrt(length(Bootstrapped_XL_Pred_Comp));
        SEM_Y = std(Bootstrapped_YL_Comp)/sqrt(length(Bootstrapped_YL_Comp));
        % Calculate the t-score for a 95% confidence interval.
        ts_X = tinv([0.025 0.975], length(Bootstrapped_XL_Pred_Comp)-1);
        ts_Y = tinv([0.025 0.975], length(Bootstrapped_YL_Comp)-1);
        % And calculate the confidence interval!
        CI_X = mean(Bootstrapped_XL_Pred_Comp) + ts_X*SEM_X;
        CI_Y = mean(Bootstrapped_YL_Comp) + ts_Y*SEM_Y;
        % Assign to the output variable
        Cognition_Onto_Eta_Gamma_PGS_PLS_CI_X(pred,comp,:) = CI_X;
        Cognition_Onto_Eta_Gamma_PGS_PLS_CI_Y(comp,:) = CI_Y;
    end
end
   
% Now, calculate the relative importance of each predictor to the PLS
% overall, regardless of components, through the VIP scores.
% Calculate the normalised beta weights from the observed PLS.
W0 = Stats_Observed.W ./ sqrt(sum(Stats_Observed.W.^2,1));
% Calculate the VIP scores for the two components.
p = size(XL_Observed,1);
sumSq = sum(XS_Observed.^2,1).*sum(YL_Observed.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
% And find variables with a VIP score greater than 1. 
indVIP = find(vipScore >= 1);
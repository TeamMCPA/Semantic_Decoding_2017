
%%% 
% Demonstration of fNIRS decoding for 8 stimulus classes
% Neurophotonics special issue featuring SfNIRS 2016 submissions
% Benjamin Zinszer, Laurie Bayet, Lauren Emberson, 
% Rajeev D.S. Raizada, & Richard Aslin
% 27 March 2017
%%%

%% Find files and define a few parameters
clear all;

DecodingResults.Between_Subjs = NaN;
DecodingResults.Semantic = NaN;

% Check in the directories Probe1 and Probe2 for *.nirs files of
% participants to include in the analysis. The search will be guided by the
% contents of the Probe1 directory.
incl_subjects = [1:24];
Probe1List = dir('Probe1/*.nirs');
Probe1List = Probe1List(incl_subjects);
Probe2List = dir('Probe2/*.nirs');
Probe2List = Probe2List(incl_subjects);

% The eight stimulus classes listed in order of appearance in the *.nirs
% file (excluding the fireworks).
conditions = {'kitty' 'bunny' 'dog' 'bear' 'foot' 'hand' 'mouth' 'nose'};

% A visual inspection of the data suggested that the most information for
% decoding existed in a window around 6.5 to 9 s after stimulus
% presentation. This corresponds to the rise of the hemodynamic response
% and the end of the trial.
scan_window = [65 90]; % In scans (collected at 10 Hz)

% Initialize a struct for storing each subject's data.
subj_struct = struct('ID',num2cell(1:length(Probe1List)));


%% Define channel stability
% Voxel stability is a typical pre-decoding analysis in functional MRI to
% identify voxels that behave similarly across runs, with no a priori
% assumptions about how voxels should respond to given stimuli (and thus no
% double-dipping). Here, we implement channel stability to the same ends.

% Initialize an array to store channel stability data
chan_selection = nan(46,length(subj_struct),2);

% Set the cutoff for stable channels. If this value is between -1 and 1, it
% will be a minimum threshold for channel stability. If this value is
% between 1 and 100, it will be a percentile cutoff (e.g., top 50th
% percentile of channels) for the given subject.
chan_stable_cut = 1; % value of 1 effectively turns off stability


%% Extract the NIRS data necessary for analysis

% Loops through each subject in the list of discovered files
for subjnum = 1:length(Probe1List),
    
    % Define files to look for
    fprintf('Subject %s:\n',Probe1List(subjnum).name);
    p1_file_path = ['Probe1/' Probe1List(subjnum).name];
    p2_file_path = ['Probe2/' Probe2List(subjnum).name];
    
    
    % This calls a function "extract_vecs_1subj" which parses the *.nirs
    % file, searches it for the various stimulus onsets, clips out the
    % epoched time series data, and saves it in a structure.
    subj_struct(subjnum).Condition = extract_vecs_1subj({p1_file_path}, conditions, scan_window, 'obsv');
    fprintf('Extracted subject %g data\n',subjnum);
    
    % Take that condition data average it down into an array with the
    % dimensions Channel X Condition
    for condnum = 1:length(subj_struct(subjnum).Condition),
        subj_struct(subjnum).rsavecs(:,condnum) = subj_struct(subjnum).Condition(condnum).ts_average;
    end
    fprintf('Saved response vectors in one array\n\n');
    
    
    % Use channel stability to select channels for inclusion. This calls
    % the function channel_stability_calculator which wraps another
    % function called StatMap_ST designed for voxel stability analysis. A
    % list of "stable channels" (defined by the cutoff value from line 44)
    % is returned, along with the stability estimates for each channel.
    
    subj_struct(subjnum).all_chans = 1:size(subj_struct(subjnum).rsavecs,1);
    
    [subj_struct(subjnum).keep_chans, stabilities]= channel_stability_calculator(subj_struct(subjnum),chan_stable_cut);
    
    chan_selection(1:length(stabilities),subjnum,1) = stabilities;
    chan_selection(subj_struct(subjnum).keep_chans,subjnum,2) = stabilities(subj_struct(subjnum).keep_chans);
    
    subj_struct(subjnum).stability = stabilities;
    
end

mean_trials_per_cond = arrayfun( @(y) ...
    mean(arrayfun(@(x) sum(subj_struct(y).Condition(x).good_trials),1:length(subj_struct(y).Condition))),...
    1:length(incl_subjects));
fprintf('Mean non-artifact trials per condition in each subject:\n');
fprintf('%0.1f ',mean_trials_per_cond);
fprintf('\n');


%% Perform between-subjects decoding using RSA
% Using only the stable channels (keep_chans) for each subject, build the
% RSA matrices for each subject and save them in a 3-D matrix of dimensions
% NumClasses X NumClasses X NumSubjects
sim_struct_indv = nan(length(conditions),length(conditions),length(subj_struct));
for i = 1:length(subj_struct),
    sim_struct_indv(:,:,i) = atanh(corr(subj_struct(i).rsavecs(subj_struct(i).keep_chans,:),'rows','pairwise'));
end

 sim_struct(:,:,1) = mean(sim_struct_indv(:,:,1:8),3);
 sim_struct(:,:,2) = mean(sim_struct_indv(:,:,9:16),3);
 sim_struct(:,:,3) = mean(sim_struct_indv(:,:,17:24),3);

% Run all pairwise comparisons between the 8 classes for each subject,
% decoding based on the group model of all remaining subjects.
% Method for this pairwise comparison is described in:
% (Anderson, Zinszer, & Raizada, 2015, Neuroimage) and on the SfNIRS poster

% This calls the function pairwise_rsa_leaveoneout which iterates through
% the subjects in an n-fold cross-validation and performs all pairwise
% comparisons for each participant (see pairwise_rsa_test)
DecodingResults.Between_Subjs = pairwise_rsa_leaveoneout(sim_struct);

% This is a quick-and-dirty parametric test of significance, which does NOT
% properly apply to cross-validation, but is convenient for the moment. We
% will run the randomization test to get a real significance test later.
[~, p , ~, stat] = ttest(DecodingResults.Between_Subjs, 0.5);
fprintf('----------------------\nQUICK AND DIRTY T-test:\n');
fprintf('Mean: %0.2f, T(%g)=%1.1f, p=%0.2f\n',nanmean(DecodingResults.Between_Subjs),stat.df,stat.tstat,p);
fprintf('DO NOT USE FOR SIGNIFICANCE TESTING\n----------------------\n\n');

% Figure 3 plots each subject's mean pairwise decoding accuracy against the
% group (other n-1 subjects)
figure
bar(DecodingResults.Between_Subjs)
title(sprintf('Between-subjects decoding accuracy, Mean: %0.2f',nanmean(DecodingResults.Between_Subjs)));
xlabel('Subject ID');


%% Load and transform the external models
% The semantic model is a 400 dimensional vector for each word from Baroni
% et al's (2014) COMPOSES corpus model. This model is state-of-the-art in
% corpus-based semantic representation and has proven powerful for fMRI
% decoding.
load('model_and_rsa_data/COMPOSES_semantic_vectors_01June2016','semantic_matrix');
semantic_sim_struct = atanh(corr(semantic_matrix));


%% Decoding each subject based on external model
% Here the external RSA model for semantic representation is applied to
% attempt decoding of the fNIRS RSA representations. Only the stable
% channels identified in the foregoing section are used here.

for subjnum = 1:size(sim_struct,3),
    
    tic;   
    fprintf('Decoding subject %g... ',subjnum);
    neural_sim_struct = sim_struct(:,:,subjnum);
    DecodingResults.Semantic(subjnum) = nanmean(pairwise_rsa_test(neural_sim_struct,semantic_sim_struct));    
    toc;
    
end

%% Plot semantic decoding
% Plots each subject's mean pairwise decoding accuracy against the
% semantic model (COMPOSES)
figure
bar( DecodingResults.Semantic );
title(sprintf('Semantic decoding accuracy, Mean: %0.2f',...
    nanmean(DecodingResults.Semantic)));
ax = gca;
xlabel('Subject ID');
ax.XTickLabel = {'Group 1', 'Group 2', 'Group3'};
fprintf('Semantic Mean: %0.2f\n',mean(DecodingResults.Semantic))

%% Group vs. group decoding
fprintf('Group vs. Group Decoding Accuracy\n');
fprintf('\tGrp1\tGrp2\tGrp3\n');
for g1 = 1:size(sim_struct,3)
    fprintf('Grp%1.0f\t',g1);
    for g2 = 1:size(sim_struct,3)
        if g2>=g1
            fprintf('%0.2f\t',mean(pairwise_rsa_test(sim_struct(:,:,g1),sim_struct(:,:,g2))));
        else
            fprintf('\t');
        end
    end
    fprintf('\n');
end

%% Null Distributions for significance testing

fprintf('\n\nBegin generating null distributions?\n');
fprintf('This procedure may take several minutes to a few hours.\n');
input('\nPress ENTER to continue... ')
fprintf('OK! Here we go!\n');

allsubj_between_subj = pairwise_rsa_leaveoneout_null(sim_struct,10000);
p_btwsubj = mean(mean(allsubj_between_subj,2)>=mean(DecodingResults.Between_Subjs));
disp(['Mean Between-Subject Decoding: ' num2str(mean(DecodingResults.Between_Subjs)) ' p=' num2str(p_btwsubj)]);

allsubj_comp_model = pairwise_rsa_null_severalsubj(sim_struct,semantic_sim_struct,10000);
p_comp = mean(mean(allsubj_comp_model,2)>=mean(DecodingResults.Semantic));
disp(['Mean COMPOSES Decoding: ' num2str(mean(DecodingResults.Semantic)) ' p=' num2str(p_comp)]);

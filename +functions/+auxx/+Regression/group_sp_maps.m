clc;
clear all;

% Import required helper functions
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*
import functions.auxx.OptimizedOperations.*

% Define parameters for analysis
prc = 90;             % Percentile threshold (e.g., for visualizing distributions)
cross_index = 0;      % Cross-validation index (0 = not using CV)
num_groups = 5;       % Number of age groups (for age stratification)
age_min = 0;          % Minimum age for inclusion
age_max = 100;        % Maximum age for inclusion
% Path to the JSON file with model result metadata
json_path = 'E:/test262025_delete/XIALPHANET.json';

% Automatically determine base directory from JSON file path
[dataset_dir, ~, ~] = fileparts(json_path);

% Load and decode dataset JSON
dataset = jsondecode(fileread(json_path));

% Set the location field automatically based on JSON file directory
dataset.Location = dataset_dir;

% Load structural model parameters
parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));
% Set up simulation parameters
Ne = parameters.Dimensions.Ne;  % Number of electrodes
Nr = parameters.Dimensions.Nr;  % Number of ROIs
Nv = parameters.Dimensions.Nr;  % Set voxel to ROI
Nw = parameters.Dimensions.Nw;  % Number of frequency bins
%

n = length(dataset.Participants);
All_Data_tmp = cell(1, n);  % combine both x and age
Mod_Matrix = cell(1, n);
ages_tmp = nan(1, n);       % to filter valid ages later

parfor i = 1:n
    participant = dataset.Participants(i);
    participant_age = participant.Age;

    if isequal(participant.Status, 'Completed') && ...
            age_min <= participant_age && participant_age <= age_max

        try
            % Load JSON metadata
            Part_Info = jsondecode(fileread(fullfile(dataset.Location, participant.SubID, participant.FileInfo)));

            % Load required matrices
            Mod_Matrix{i} = load(fullfile(dataset.Location, participant.SubID, Part_Info.Mod_Weights));

            % Alpha component
            alpha_process = load(fullfile(dataset.Location, participant.SubID, Part_Info.Alpha_estimate));
            a = zeros(length(alpha_process.Power), 4);
            a(:,1) = alpha_process.Power;
            a(:,2) = alpha_process.Width;
            a(:,3) = alpha_process.Exponent;
            a(:,4) = alpha_process.PAF;
            A(i,:,:) = a;

            % Xi component
            xi_process = load(fullfile(dataset.Location, participant.SubID, Part_Info.Xi_estimate));
            e = zeros(length(xi_process.Power), 3);
            e(:,1) = xi_process.Power;
            e(:,2) = xi_process.Width;
            e(:,3) = xi_process.Exponent;
            E(i,:,:) = e;

            % Save the age
            ages_tmp(i) = participant_age;

        catch ME
            warning("Participant %s failed: %s", participant.SubID, ME.message);
        end
    end
end

%% === Compute group averages for Xi (E) and Alpha (A) ===

% Define age bins (0–100, 5 groups of 20 years each)
edges = linspace(0, 100, 6);   % [0 20 40 60 80 100]
num_groups = length(edges) - 1;

% Make sure output directory exists
outdir = 'E:/group_average_xialphanet';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Initialize containers for all groups
Xi_Power_all     = cell(1,num_groups);
Xi_Bandwidth_all = cell(1,num_groups);
Xi_Exponent_all  = cell(1,num_groups);

Alpha_Power_all     = cell(1,num_groups);
Alpha_Bandwidth_all = cell(1,num_groups);
Alpha_Exponent_all  = cell(1,num_groups);
Alpha_Frequency_all = cell(1,num_groups);

for g = 1:num_groups
    % Define age range for group g
    age_min_g = edges(g);
    age_max_g = edges(g+1);

    % Select participants in this age bin
    idx = ages_tmp >= age_min_g & ages_tmp <= age_max_g;

    % Compute group means if there are valid participants
    if any(idx)
        % Xi parameters
        Xi_Power_all{g}     = squeeze(mean(E(idx,:,1), 1, 'omitnan')); 
        Xi_Bandwidth_all{g} = squeeze(mean(E(idx,:,2), 1, 'omitnan')); 
        Xi_Exponent_all{g}  = squeeze(mean(E(idx,:,3), 1, 'omitnan')); 

        % Alpha parameters
        Alpha_Power_all{g}     = squeeze(mean(A(idx,:,1), 1, 'omitnan'));
        Alpha_Bandwidth_all{g} = squeeze(mean(A(idx,:,2), 1, 'omitnan'));
        Alpha_Exponent_all{g}  = squeeze(mean(A(idx,:,3), 1, 'omitnan'));
        Alpha_Frequency_all{g} = squeeze(mean(A(idx,:,4), 1, 'omitnan'));
    else
        warning('No participants in age group %d (%d–%d).', g, age_min_g, age_max_g);
        Xi_Power_all{g}     = nan(size(E,2),1);
        Xi_Bandwidth_all{g} = nan(size(E,2),1);
        Xi_Exponent_all{g}  = nan(size(E,2),1);
        Alpha_Power_all{g}     = nan(size(A,2),1);
        Alpha_Bandwidth_all{g} = nan(size(A,2),1);
        Alpha_Exponent_all{g}  = nan(size(A,2),1);
        Alpha_Frequency_all{g} = nan(size(A,2),1);
    end
end

%% === Save per parameter ===

% Xi Power
for g = 1:num_groups
    eval(sprintf('Xi_Power_group%d = Xi_Power_all{g};', g));
end
save(fullfile(outdir,'Xi_Power.mat'), 'Xi_Power_group1','Xi_Power_group2','Xi_Power_group3','Xi_Power_group4','Xi_Power_group5');

% Xi Bandwidth
for g = 1:num_groups
    eval(sprintf('Xi_Bandwidth_group%d = Xi_Bandwidth_all{g};', g));
end
save(fullfile(outdir,'Xi_Bandwidth.mat'), 'Xi_Bandwidth_group1','Xi_Bandwidth_group2','Xi_Bandwidth_group3','Xi_Bandwidth_group4','Xi_Bandwidth_group5');

% Xi Exponent
for g = 1:num_groups
    eval(sprintf('Xi_Exponent_group%d = Xi_Exponent_all{g};', g));
end
save(fullfile(outdir,'Xi_Exponent.mat'), 'Xi_Exponent_group1','Xi_Exponent_group2','Xi_Exponent_group3','Xi_Exponent_group4','Xi_Exponent_group5');

% Alpha Power
for g = 1:num_groups
    eval(sprintf('Alpha_Power_group%d = Alpha_Power_all{g};', g));
end
save(fullfile(outdir,'Alpha_Power.mat'), 'Alpha_Power_group1','Alpha_Power_group2','Alpha_Power_group3','Alpha_Power_group4','Alpha_Power_group5');

% Alpha Bandwidth
for g = 1:num_groups
    eval(sprintf('Alpha_Bandwidth_group%d = Alpha_Bandwidth_all{g};', g));
end
save(fullfile(outdir,'Alpha_Bandwidth.mat'), 'Alpha_Bandwidth_group1','Alpha_Bandwidth_group2','Alpha_Bandwidth_group3','Alpha_Bandwidth_group4','Alpha_Bandwidth_group5');

% Alpha Exponent
for g = 1:num_groups
    eval(sprintf('Alpha_Exponent_group%d = Alpha_Exponent_all{g};', g));
end
save(fullfile(outdir,'Alpha_Exponent.mat'), 'Alpha_Exponent_group1','Alpha_Exponent_group2','Alpha_Exponent_group3','Alpha_Exponent_group4','Alpha_Exponent_group5');

% Alpha Frequency
for g = 1:num_groups
    eval(sprintf('Alpha_Frequency_group%d = Alpha_Frequency_all{g};', g));
end
save(fullfile(outdir,'Alpha_Frequency.mat'), 'Alpha_Frequency_group1','Alpha_Frequency_group2','Alpha_Frequency_group3','Alpha_Frequency_group4','Alpha_Frequency_group5');

disp('Group averages saved per parameter.');



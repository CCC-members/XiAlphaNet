function compute_group_averages(json_path, outdir, structural_dir, num_groups, age_min, age_max)
% COMPUTE_GROUP_AVERAGES
% -------------------------------------------------------------------------
% This function computes group-averaged Xi and Alpha spectral parameters 
% from participant-level Xi-AlphaNET model outputs. 
%
% INPUTS:
%   json_path     : Path to the dataset JSON file (metadata of participants).
%   outdir        : Output directory where results will be saved.
%   structural_dir: Directory containing the structural parameters.mat file.
%   num_groups    : Number of age groups for binning participants.
%   age_min       : Minimum age for inclusion (default: 0).
%   age_max       : Maximum age for inclusion (default: 100).
%
% OUTPUT:
%   Saves .mat files with group-averaged parameter maps per age bin:
%   - Xi_Power.mat, Xi_Bandwidth.mat, Xi_Exponent.mat
%   - Alpha_Power.mat, Alpha_Bandwidth.mat, Alpha_Exponent.mat, Alpha_Frequency.mat
%
% Example usage:
%   compute_group_averages('E:/XIALPHANET.json', ...
%                          'E:/group_average_xialphanet', ...
%                          'E:/structural', 5, 0, 100);
% -------------------------------------------------------------------------

%% === SETUP ===
if nargin < 4, num_groups = 5; end
if nargin < 5, age_min = 0; end
if nargin < 6, age_max = 100; end

% Ensure output directory exists
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Load dataset metadata
dataset_dir = fileparts(json_path);
dataset = jsondecode(fileread(json_path));
dataset.Location = dataset_dir;

% Load structural parameters
parameters = load(fullfile(structural_dir, 'parameters.mat'));
Ne = parameters.Dimensions.Ne; %#ok<NASGU>
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nr; %#ok<NASGU>
Nw = parameters.Dimensions.Nw; %#ok<NASGU>

%% === INITIALIZE ===
n = length(dataset.Participants);
A = []; % Alpha parameters
E = []; % Xi parameters
ages_tmp = nan(1,n);

% Preallocate containers
All_Data_tmp = cell(1, n); %#ok<NASGU>
Mod_Matrix   = cell(1, n);

%% === LOOP OVER PARTICIPANTS ===
parfor i = 1:n
    participant = dataset.Participants(i);
    participant_age = participant.Age;

    if isequal(participant.Status, 'Completed') && ...
            age_min <= participant_age && participant_age <= age_max

        try
            % Load participant-specific metadata
            Part_Info = jsondecode(fileread(fullfile(dataset.Location, ...
                                                     participant.SubID, ...
                                                     participant.FileInfo)));

            % Load weights
            Mod_Matrix{i} = load(fullfile(dataset.Location, participant.SubID, Part_Info.Mod_Weights));

            % Alpha parameters: Power, Bandwidth, Exponent, Frequency
            alpha_process = load(fullfile(dataset.Location, participant.SubID, Part_Info.Alpha_estimate));
            a = [alpha_process.Power, alpha_process.Width, ...
                 alpha_process.Exponent, alpha_process.PAF];
            A(i,:,:) = a;

            % Xi parameters: Power, Bandwidth, Exponent
            xi_process = load(fullfile(dataset.Location, participant.SubID, Part_Info.Xi_estimate));
            e = [xi_process.Power, xi_process.Width, xi_process.Exponent];
            E(i,:,:) = e;

            % Save participant age
            ages_tmp(i) = participant_age;

        catch ME
            warning("Participant %s failed: %s", participant.SubID, ME.message);
        end
    end
end

%% === DEFINE AGE BINS ===
edges = linspace(age_min, age_max, num_groups+1);

% Containers for results
Xi_Power_all     = cell(1,num_groups);
Xi_Bandwidth_all = cell(1,num_groups);
Xi_Exponent_all  = cell(1,num_groups);

Alpha_Power_all     = cell(1,num_groups);
Alpha_Bandwidth_all = cell(1,num_groups);
Alpha_Exponent_all  = cell(1,num_groups);
Alpha_Frequency_all = cell(1,num_groups);

%% === GROUP AVERAGING ===
for g = 1:num_groups
    % Select participants in this age range
    idx = ages_tmp >= edges(g) & ages_tmp < edges(g+1);

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
        warning('No participants in age group %d (%d–%d).', g, edges(g), edges(g+1));
        Xi_Power_all{g}        = nan(Nr,1);
        Xi_Bandwidth_all{g}    = nan(Nr,1);
        Xi_Exponent_all{g}     = nan(Nr,1);
        Alpha_Power_all{g}     = nan(Nr,1);
        Alpha_Bandwidth_all{g} = nan(Nr,1);
        Alpha_Exponent_all{g}  = nan(Nr,1);
        Alpha_Frequency_all{g} = nan(Nr,1);
    end
end

%% === SAVE RESULTS ===
save(fullfile(outdir,'Xi_Power.mat'), 'Xi_Power_all');
save(fullfile(outdir,'Xi_Bandwidth.mat'), 'Xi_Bandwidth_all');
save(fullfile(outdir,'Xi_Exponent.mat'), 'Xi_Exponent_all');
save(fullfile(outdir,'Alpha_Power.mat'), 'Alpha_Power_all');
save(fullfile(outdir,'Alpha_Bandwidth.mat'), 'Alpha_Bandwidth_all');
save(fullfile(outdir,'Alpha_Exponent.mat'), 'Alpha_Exponent_all');
save(fullfile(outdir,'Alpha_Frequency.mat'), 'Alpha_Frequency_all');

fprintf('Group averages successfully saved in %s\n', outdir);

end

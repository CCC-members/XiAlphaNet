%% Population spectral net
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
num_groups = 1;       % Number of age groups (for age stratification)
mode = 1;             % 0 = amplitude plots, 1 = zero-inflation probability plots
age_min = 0;          % Minimum age for inclusion
age_max = 100;        % Maximum age for inclusion
Nw = 47;
% Path to the JSON file with model result metadata. Modify this directions
% manually acording to the location of the downloaded data
json_path = 'D:\OneDrive - CCLAB\New_Data_XiAlphaNET\xialphanet_newresults22\XIALPHANET.json';
dir_data = 'D:\OneDrive - CCLAB\New_Data_XiAlphaNET\HarMNqEEG_norms\MultinationalNorms';

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

%%
%% Load template cortex and ROI projection matrix
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit

Cortex   = load("templates/Cortex.mat");
template = load("templates/axes.mat");
currentAxes = template.axes;
hemis = [];
[Cortex, iHideVert] = split_hemisphere(Cortex, hemis);

% ROI projection matrix (Nv × 360)
[R,R_inv] = functions.auxx.ModelVectorization.roi_average_operators(Cortex,10);
%%

n = length(dataset.Participants);
All_Data_tmp = cell(1, n);  % combine both x and age
Mod_Matrix = cell(1, n);
ages_tmp = nan(1, n);       % to filter valid ages later



parfor i = 1:n
    i
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
            ta = set_threshold_em(a(:,1));
            a(:,1) = a(:,1) .* (a(:,1) > ta);
            a(:,2) = alpha_process.Width;
            a(:,3) = alpha_process.Exponent;
            a(:,4) = alpha_process.PAF;
            ta = set_threshold_em(a(:,1));
            a(:,1) = a(:,1) .* (a(:,1) > ta);

            % Xi component
            xi_process = load(fullfile(dataset.Location, participant.SubID, Part_Info.Xi_estimate));
            e = zeros(length(xi_process.Power), 3);
            e(:,1) = xi_process.Power;
            e(:,2) = xi_process.Width;
            e(:,3) = xi_process.Exponent;
            te = set_threshold_em(e(:,1));
            e(:,1) = e(:,1) .* (e(:,1) > te);

            % Map voxel to roi
            eR = zeros(size(R,1),3);
            aR = zeros(size(R,1),4);

            eR(:,1)  = R * e(:,1);
            eR(:,2)  = R * e(:,2);
            eR(:,3)  = R * e(:,3);

            aR(:,1)  = R * a(:,1);
            aR(:,2)  = R * a(:,2);
            aR(:,3)  = R * a(:,3);
            aR(:,4)  = R * a(:,4);

            % Build the combined feature
            x = v2x(eR, aR, 0.05);

            % Assign single cell entry
            All_Data_tmp{i} = {x, participant_age};
            ages_tmp(i) = participant_age;

        catch ME
            warning("Participant %s failed: %s", participant.SubID, ME.message);
        end
    end
end

% Post-processing to unpack All_Data
valid = ~cellfun(@isempty, All_Data_tmp);
All_Data = cell(2, sum(valid));
All_Data(1,:) = cellfun(@(c) c{1}, All_Data_tmp(valid), 'UniformOutput', false);
All_Data(2,:) = cellfun(@(c) c{2}, All_Data_tmp(valid), 'UniformOutput', false);

Delay_Matrix = Mod_Matrix(valid);
ages = ages_tmp(valid);

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
% Path to the JSON file with model result metadata
json_path = '/mnt/Store/Ronaldo/dev/Data/NewFolder/XIALPHANET.json';

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

            % Build the combined feature
            x = v2x(e, a, 1);

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

%%

% Set data directory for simulation
dir_data = '/mnt/Store/Ronaldo/dev/Data/norms';
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
selected_folders = subject_folders(randperm(length(subject_folders), 1));
subject_folder = selected_folders(1).name;
mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
data_struct = load(mat_file_path);
freq = data_struct.data_struct.freqrange(1:Nw);
parameters.Data.freq = freq;


% Set model parameters for Xi-AlphaNET estimation
properties.model_params.nFreqs = Nw;
properties.general_params.parallel.conn_delay = 1;

%% === Setup ===
R = parameters.Model.R;
C0 = parameters.Model.C;
D0 = parameters.Model.D;
Nr = size(R, 1);
N = length(ages);  % Total number of participants

batchSize = 20;
nBatches = ceil(N / batchSize);

% --- Initialize overall sums ---
c_xi = zeros(Nr);
c_alpha = zeros(Nr);

% --- Define 5 age groups: 020, 2040, 4060, 6080, 80100 ---
ageEdges = [0, 20, 40, 60, 80, 101];
nGroups = length(ageEdges) - 1;
ageGroupLabels = {'020', '2040', '4060', '6080', '80100'};

c_xi_group = cell(1, nGroups);
c_alpha_group = cell(1, nGroups);
count_group = zeros(1, nGroups);

for gg = 1:nGroups
    c_xi_group{gg} = zeros(Nr);
    c_alpha_group{gg} = zeros(Nr);
end

% === Main Batch Loop ===
for b = 1:nBatches
    b
    batchStart = (b - 1) * batchSize + 1;
    batchEnd = min(b * batchSize, N);
    idx = batchStart:batchEnd;
    nThisBatch = length(idx);

    result(nThisBatch) = struct('cxi', zeros(Nr), 'calpha', zeros(Nr), 'group', 0);

    parfor jj = 1:nThisBatch
        j = idx(jj)
        x = All_Data{1, j};
        [e, a, s2] = x2v(x);

        omega_xi = freq(1);
        omega_alpha = freq(25);

        C = Mod_Matrix{j}.Mod_Weights(2) * C0;
        D = Mod_Matrix{j}.Mod_Weights(1) * D0;
        I = zeros(size(C));

        Z_xi    = I + C .* exp(-2 * pi * 1i * omega_xi * D);
        Z_alpha = I + C .* exp(-2 * pi * 1i * omega_alpha * D);

        Tj_cross_xi    = R * Z_xi;
        Tj_cross_alpha = R * Z_alpha;

        xi_omega = e(:,1) ./ (1 + e(:,2) .* omega_xi.^2).^e(:,3);
        alpha_omega = 0.5 * (a(:,1) ./ (1 + a(:,2) .* (omega_alpha - a(:,4)).^2).^a(:,3) + ...
                             a(:,1) ./ (1 + a(:,2) .* (omega_alpha + a(:,4)).^2).^a(:,3));

        % Compute connectivity and log-transform
        result(jj).cxi = logm(computeTDT(Tj_cross_xi, xi_omega) + 1e-6 * eye(Nr));
        result(jj).calpha = logm(computeTDT(Tj_cross_alpha, alpha_omega) + 1e-6 * eye(Nr));

        % Assign age group
        age_j = ages(j);
        group_j = find((age_j >= ageEdges(1:end-1)) & (age_j < ageEdges(2:end)), 1);
        result(jj).group = group_j;
    end

    % Accumulate overall and group-wise results
    for jj = 1:nThisBatch
        gg = result(jj).group;
        if ~isempty(gg)
            c_xi = c_xi + result(jj).cxi;
            c_alpha = c_alpha + result(jj).calpha;

            c_xi_group{gg} = c_xi_group{gg} + result(jj).cxi;
            c_alpha_group{gg} = c_alpha_group{gg} + result(jj).calpha;
            count_group(gg) = count_group(gg) + 1;
        end
    end
end

% === Finalize Averages ===
c_xi = expm(c_xi / N);
C_alpha = expm(c_alpha / N);

C_xi_groups = cell(1, nGroups);
C_alpha_groups = cell(1, nGroups);
for gg = 1:nGroups
    if count_group(gg) > 0
        C_xi_groups{gg} = expm(c_xi_group{gg} / count_group(gg));
        C_alpha_groups{gg} = expm(c_alpha_group{gg} / count_group(gg));
    else
        C_xi_groups{gg} = nan(Nr);
        C_alpha_groups{gg} = nan(Nr);
    end
end

%% === Visualization ===
% Overall Xi Process
A = c_xi;
guide.Visualization.effective_spectral_net;
title('Overall Effective Connectivity: Xi Process');

% Overall Alpha Process
A = C_alpha;
guide.Visualization.effective_spectral_net;
title('Overall Effective Connectivity: Alpha Process');

% Age Group Visualizations
for gg = 1:nGroups
    gg
    ageLabel = ['Age Group ' ageGroupLabels{gg}];

    if ~all(isnan(C_xi_groups{gg}(:)))
        A = C_xi_groups{gg};
        guide.Visualization.effective_spectral_net;
       % title(['Xi Process - ' ageLabel]);
    end

    if ~all(isnan(C_alpha_groups{gg}(:)))
        A = C_alpha_groups{gg};
        guide.Visualization.effective_spectral_net;
        %title(['Alpha Process - ' ageLabel]);
    end
end

guide.Visualization.effective_spectral_net

%%


clc;
clear all;

% Import required helper functions
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*

% Define parameters for analysis
prc = 90;             % Percentile threshold (e.g., for visualizing distributions)
cross_index = 0;      % Cross-validation index (0 = not using CV)
num_groups = 1;       % Number of age groups (for age stratification)
mode = 0;             % 0 = amplitude plots, 1 = zero-inflation probability plots
age_min = 0;          % Minimum age for inclusion
age_max = 100;        % Maximum age for inclusion

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
Cortex                  = load("templates/Cortex.mat");




%% Positive colormap: white ? orange ? dark orange
n_gray = 95;     % Gray range from 0 to 0.95
n_color = 5;     % Color range from 0.95 to 1
gray = [0.85, 0.85, 0.85];

clip01 = @(x) min(max(x, 0), 1);

% Colors
dark_orange   = clip01([0.60, 0.30, 0.00]);
bright_orange = clip01([1.00, 0.55, 0.00]);

dark_green    = clip01([0.00, 0.35, 0.10]);
bright_green  = clip01([0.10, 0.85, 0.50]);

dark_blue     = clip01([0.15, 0.30, 0.65]);
bright_blue   = clip01([0.30, 0.60, 1.00]);

% Create gray part (flat gray)
gray_part = repmat(gray, n_gray, 1);

% Positive transitions
orange_positive = [gray_part; ...
                   [linspace(gray(1), bright_orange(1), n_color)', ...
                    linspace(gray(2), bright_orange(2), n_color)', ...
                    linspace(gray(3), bright_orange(3), n_color)']];

green_positive = [gray_part; ...
                  [linspace(gray(1), bright_green(1), n_color)', ...
                   linspace(gray(2), bright_green(2), n_color)', ...
                   linspace(gray(3), bright_green(3), n_color)']];

blue_positive = [gray_part; ...
                 [linspace(gray(1), bright_blue(1), n_color)', ...
                  linspace(gray(2), bright_blue(2), n_color)', ...
                  linspace(gray(3), bright_blue(3), n_color)']];

% Negative transitions (also gray to dark color)
orange_negative = [gray_part; ...
                   [linspace(gray(1), dark_orange(1), n_color)', ...
                    linspace(gray(2), dark_orange(2), n_color)', ...
                    linspace(gray(3), dark_orange(3), n_color)']];

green_negative = [gray_part; ...
                  [linspace(gray(1), dark_green(1), n_color)', ...
                   linspace(gray(2), dark_green(2), n_color)', ...
                   linspace(gray(3), dark_green(3), n_color)']];

blue_negative = [gray_part; ...
                 [linspace(gray(1), dark_blue(1), n_color)', ...
                  linspace(gray(2), dark_blue(2), n_color)', ...
                  linspace(gray(3), dark_blue(3), n_color)']];

% Optional preview
figure;
subplot(3,2,1); imagesc(permute(orange_negative, [1 3 2])); title('Orange Negative');
subplot(3,2,2); imagesc(permute(orange_positive, [1 3 2])); title('Orange Positive');

subplot(3,2,3); imagesc(permute(green_negative, [1 3 2])); title('Green Negative');
subplot(3,2,4); imagesc(permute(green_positive, [1 3 2])); title('Green Positive');

subplot(3,2,5); imagesc(permute(blue_negative, [1 3 2]));  title('Blue Negative');
subplot(3,2,6); imagesc(permute(blue_positive, [1 3 2]));  title('Blue Positive');



%%

ages = [];
All_Data = {}; 
index = 1;
HCP_MMP1 = Cortex.Atlas(Cortex.iAtlas);
N_voxel = length(Cortex.Curvature);
N_roi = 360;
for i=1:length(dataset.Participants)
    participant = dataset.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed')) && age_min <= participant_age && participant_age <= age_max
       ages = [ages,participant_age];
       All_Data{2,index} =  participant_age;
       Part_Info = jsondecode(fileread(fullfile(dataset.Location,participant.SubID,participant.FileInfo)));
       alpha_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
       a(:,1) = alpha_process.Power;
       a(:,2) = alpha_process.Width;
       a(:,3) = alpha_process.Exponent;
       a(:,4) = alpha_process.PAF;
       xi_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
       e(:,1) = xi_process.Power;
       e(:,2) = xi_process.Width;
       e(:,3) = xi_process.Exponent;
       s2 = 1;
       x = v2x(e,a,s2);
       All_Data{1,index} = x;
       index = index +1;
    end
end

% Initialize storage for Peak Alpha Frequency (PAF), Amplitude of the Alpha, and Amplitude of Xi
% Initialize storage
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

threshold_PAF = 7;
import functions.auxx.Refine_Solution.*

% Temporary storage for thresholds
threshold_Alpha = zeros(1, length(All_Data));
threshold_Xi = zeros(1, length(All_Data));

% First pass: extract features
for j = 1:length(All_Data)
    fprintf('Processing subject %d\n', j);
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Alpha Amplitude
    XiAmp_all(:,j) = e(:,1);          % Xi Amplitude
end

% Compute mean across voxels for each subject
mean_PAF = mean(PAF_all, 1, 'omitnan');
mean_Alpha = mean(AlphaAmp_all, 1, 'omitnan');
mean_Xi = mean(XiAmp_all, 1, 'omitnan');

% Z-score computation
z_PAF = zscore(mean_PAF);
z_Alpha = zscore(mean_Alpha);
z_Xi = zscore(mean_Xi);

% Identify outliers: abs(z) > 3 in any metric
outlier_idx = abs(z_PAF) > 3 | abs(z_Alpha) > 3 | abs(z_Xi) > 3;
valid_idx = ~outlier_idx;
fprintf('Removed %d outlier subjects\n', sum(outlier_idx));

% Filter out outliers
PAF_all = PAF_all(:, valid_idx);
AlphaAmp_all = AlphaAmp_all(:, valid_idx);
XiAmp_all = XiAmp_all(:, valid_idx);
All_Data = All_Data(:, valid_idx);

% Second pass: Apply thresholds
parfor j = 1:length(All_Data)
    [e, a, s2] = x2v(All_Data{1,j});
    Alpha_j = a(:,1);
    Xi_j = e(:,1);
    PAF_j = a(:,4);

    threshold_Alpha(j) = set_threshold_em(Alpha_j);
    threshold_Xi(j) = set_threshold_em(Xi_j);
    if mode == 1 % For probability
        Alpha_j =  (Alpha_j > threshold_Alpha(j));
        Xi_j_bin = Xi_j > threshold_Alpha(j);
        PAF_j_bin = (PAF_j .* (Alpha_j > threshold_Alpha(j))) > threshold_PAF;
    else % For amplitude Distribution
        Alpha_j =  Alpha_j.*(Alpha_j > threshold_Alpha(j));
        Xi_j_bin = Xi_j.*(Xi_j > threshold_Alpha(j));
        PAF_j_bin = PAF_j.*((PAF_j .* (Alpha_j > threshold_Alpha(j))) > threshold_PAF);
    end

    AlphaAmp_all(:,j) = Alpha_j;
    XiAmp_all(:,j) = Xi_j_bin;
    PAF_all(:,j) = PAF_j_bin;
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age


%% === COMPONENT CONFIGURATION ===
components     = {'PAF', 'AlphaAmp', 'XiAmp'};
data_all       = {PAF_all, AlphaAmp_all, XiAmp_all};
colormaps_all  = {
    orange_positive, orange_negative; 
    blue_positive,   blue_negative; 
    green_positive,  green_negative
};

% === PROCESS EACH COMPONENT (VOXEL-WISE) ===
for comp_idx = 1:3
    fprintf('\nProcessing %s (voxel-wise)\n', components{comp_idx});
    current_data = data_all{comp_idx};
    cmap_pos     = colormaps_all{comp_idx, 1};
    cmap_neg     = colormaps_all{comp_idx, 2};

    % === Fit ZIG model for each voxel ===
    result_all = cell(N_voxel, 1);
    parfor v = 1:N_voxel
        Y_v = current_data(v, :);       % [1 x N] vector for voxel v
        result_all{v} = fitZIG_random(Y_v, ages);
    end

    % === GAMMA > 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.gamma(2) > 0
            J(v) = 1 - result_all{v}.p_gamma(2);
        end
    end
    guide.Visualization.esi_plot_single;
    colormap(cmap_pos); caxis([0 1]); title([components{comp_idx}, ' ? > 0']);

    % === GAMMA < 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.gamma(2) < 0
            J(v) = 1 - result_all{v}.p_gamma(2);
        end
    end
    guide.Visualization.esi_plot_single;
    colormap(cmap_neg); caxis([0 1]); title([components{comp_idx}, ' ? < 0']);

    % === BETA > 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.beta(3) > 0
            J(v) = 1 - result_all{v}.p_beta(3);
        end
    end
    guide.Visualization.esi_plot_single;
    colormap(cmap_pos); caxis([0 1]); title([components{comp_idx}, ' ß > 0']);

    % === BETA < 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.beta(3) < 0
            J(v) = 1 - result_all{v}.p_beta(3);
        end
    end
    guide.Visualization.esi_plot_single;
    colormap(cmap_neg); caxis([0 1]); title([components{comp_idx}, ' ß < 0']);
end

%%
figure 
hold on 
for  j=1:N_voxel
    plot(ages(:),result_all{j}.p_beta(3),'.')
end
hold of 
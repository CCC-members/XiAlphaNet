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

%%
%%
n = 256;  % Number of color steps

%% Positive colormap: white ? orange ? dark orange
n = 256;
n1 = floor(n/2);
n2 = n - n1;

% ORANGE (positive)  soft gray ? orange ? dark orange
gray_to_orange = [linspace(0.9, 1, n1)', linspace(0.9, 0.5, n1)', linspace(0.9, 0, n1)'];
orange_to_dark = [linspace(1, 0.6, n2)', linspace(0.5, 0.2, n2)', linspace(0, 0, n2)'];
orange_positive = [gray_to_orange; orange_to_dark];

% PURPLE (negative for orange)  soft gray ? purple ? dark purple
gray_to_purple = [linspace(0.9, 0.7, n1)', linspace(0.9, 0.5, n1)', linspace(0.9, 0.9, n1)'];
purple_to_dark = [linspace(0.7, 0.4, n2)', linspace(0.5, 0.0, n2)', linspace(0.9, 0.4, n2)'];
orange_negative = [gray_to_purple; purple_to_dark];

% BLUE (negative)  soft gray ? blue ? dark blue
gray_to_blue = [linspace(0.9, 0.4, n1)', linspace(0.9, 0.6, n1)', linspace(0.9, 1, n1)'];
blue_to_dark = [linspace(0.4, 0, n2)', linspace(0.6, 0, n2)', linspace(1, 0.4, n2)'];
blue_negative = [gray_to_blue; blue_to_dark];

% RED (optional positive for blue, not requested but added for balance)
gray_to_red = [linspace(0.9, 1, n1)', linspace(0.9, 0.4, n1)', linspace(0.9, 0.4, n1)'];
red_to_dark = [linspace(1, 0.6, n2)', linspace(0.4, 0.1, n2)', linspace(0.4, 0.1, n2)'];
blue_positive = [gray_to_red; red_to_dark];

% GREEN (positive)  soft gray ? teal green ? dark green
gray_to_green = [linspace(0.9, 0, n1)', linspace(0.9, 0.6, n1)', linspace(0.9, 0.5, n1)'];
green_to_dark = [linspace(0, 0, n2)', linspace(0.6, 0.3, n2)', linspace(0.5, 0.2, n2)'];
green_positive = [gray_to_green; green_to_dark];

% MAGENTA (negative for green)  soft gray ? magenta ? dark magenta
gray_to_magenta = [linspace(0.9, 0.9, n1)', linspace(0.9, 0.5, n1)', linspace(0.9, 0.9, n1)'];
magenta_to_dark = [linspace(0.9, 0.4, n2)', linspace(0.5, 0.0, n2)', linspace(0.9, 0.4, n2)'];
green_negative = [gray_to_magenta; magenta_to_dark];



%%

ages = [];
All_Data = {}; 
index = 1;
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
outlier_idx = abs(z_PAF) > 2.5 | abs(z_Alpha) > 2.5 | abs(z_Xi) > 2.5;
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
%%

HCP_MMP1 = Cortex.Atlas(Cortex.iAtlas);
N_voxel = length(Cortex.Curvature);
N_roi = 360;

%% PAF
% Initialize R with zeros
J = zeros(N_voxel, 1);
result_all = cell(N_roi,1);
parfor i = 1:360
    i
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    Y_roi = PAF_all(vertices_i,:);
    result = fitZIG_random(Y_roi, ages);
    result_all{i} = result; 
end

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_gamma(2))*((result_all{i}.gamma(2))>0);%*sign(0.05-mean(result_all{i}.pi_est));
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(positive);

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_gamma(2))*((result_all{i}.gamma(2))<0);%*sign(0.05-mean(result_all{i}.pi_est));
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(negative);


for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_beta(3))*((result_all{i}.beta(3))>0);
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(positive);

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_beta(3))*((result_all{i}.beta(3))<0);
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(negative);



% 
% for i=1:360
%     vertices_i = HCP_MMP1.Scouts(i).Vertices;
%     J(vertices_i) = mean((1- result_all{i}.pi_est).*(result_all{i}.mu_est));
% end
% %J = J.*(J>0.99999)
% guide.Visualization.esi_plot_single
% colormap(custom_colormap);




%% AlphaAMp
% Initialize R with zeros
J = zeros(N_voxel, 1);
result_all = cell(N_roi,1);
parfor i = 1:360
    i
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    Y_roi = AlphaAmp_all(vertices_i,:);
    result = fitZIG_random(Y_roi, ages);
    result_all{i} = result; 
end


for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_gamma(2))*((result_all{i}.gamma(2))>0);%*sign(0.05-mean(result_all{i}.pi_est));
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(positive);

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_gamma(2))*((result_all{i}.gamma(2))<0);%*sign(0.05-mean(result_all{i}.pi_est));
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(negative);


for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_beta(3))*((result_all{i}.beta(3))>0);
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(positive);

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_beta(3))*((result_all{i}.beta(3))<0);
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(negative);

% for i=1:360
%     vertices_i = HCP_MMP1.Scouts(i).Vertices;
%     J(vertices_i) = mean((1- result_all{i}.pi_est).*(result_all{i}.mu_est));
% end
% %J = J.*(J>0.99999)
% guide.Visualization.esi_plot_single
% colormap(custom_colormap);

%%


% Initialize R with zeros
J = zeros(N_voxel, 1);
result_all = cell(N_roi,1);
parfor i = 1:360
    i
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    Y_roi = XiAmp_all(vertices_i,:);
    result = fitZIG_random(Y_roi, ages);
    result_all{i} = result; 
end


for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_gamma(2))*((result_all{i}.gamma(2))>0);%*sign(0.05-mean(result_all{i}.pi_est));
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(positive);

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_gamma(2))*((result_all{i}.gamma(2))<0);%*sign(0.05-mean(result_all{i}.pi_est));
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(negative);


for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_beta(3))*((result_all{i}.beta(3))>0);
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(positive);

for i=1:360
    vertices_i = HCP_MMP1.Scouts(i).Vertices;
    J(vertices_i) = (1- result_all{i}.p_beta(3))*((result_all{i}.beta(3))<0);
end
%J = J.*(J>0.99999)
guide.Visualization.esi_plot_single
colormap(negative);


%%
% === Setup ===
signal_types = {'PAF', 'AlphaAmp', 'XiAmp'};
signal_data = {PAF_all, AlphaAmp_all, XiAmp_all};
N_roi = length(HCP_MMP1.Scouts);

for s = 1:numel(signal_types)
    signal_name = signal_types{s};
    fprintf('--- Processing %s ---\n', signal_name);
    current_signal = signal_data{s};

    % Initialize outputs
    J = zeros(N_voxel, 1);
    result_all = cell(N_roi, 1);

    % === ZIG model fitting across ROIs ===
    parfor i = 1:N_roi
        fprintf('Fitting ZIG: ROI %d\n', i);
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        Y_roi = current_signal(vertices_i, :);
        result_all{i} = fitZIG_random(Y_roi, ages);
    end

    % === Visualization 1: Pi Significance (p_gamma) ===
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        pi_mean = mean(result_all{i}.pi_est);
        J(vertices_i) = (1 - result_all{i}.p_gamma(2)) * sign(0.5 - pi_mean);
    end
    guide.Visualization.esi_plot_single;
    title(sprintf('%s  Significance of Zero-Inflation (p\\_\\gamma)', signal_name), 'Interpreter', 'none');
    colormap(custom_colormap);

    % === Visualization 2: Beta Significance (p_beta) ===
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        J(vertices_i) = (1 - result_all{i}.p_beta(3)) * sign(result_all{i}.beta(3));
    end
    guide.Visualization.esi_plot_single;
    title(sprintf('%s  Significance of Age Slope (p\\_\\beta)', signal_name), 'Interpreter', 'none');
    colormap(custom_colormap);

    % === Visualization 3: Weighted Mean of Signal ((1 - pi) * mu) ===
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        J(vertices_i) = mean((1 - result_all{i}.pi_est) .* result_all{i}.mu_est);
    end
    guide.Visualization.esi_plot_single;
    title(sprintf('%s  Weighted Mean Activity ((1 - \\pi) * \\mu)', signal_name), 'Interpreter', 'tex');
    colormap(custom_colormap);
end

%%
figure;

var_names = {'PAF', 'AlphaAmp', 'XiAmp'};
data_all = {PAF_all, AlphaAmp_all, XiAmp_all};
pos_cmaps = {orange_positive, blue_positive, green_positive};
neg_cmaps = {orange_negative, blue_negative, green_negative};

for var_idx = 1:3
    % Initialize
    J = zeros(N_voxel, 1);
    result_all = cell(N_roi,1);
    current_data = data_all{var_idx};

    % Run ZIG per ROI
    parfor i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        Y_roi = current_data(vertices_i,:);
        result_all{i} = fitZIG_random(Y_roi, ages); 
    end

    % --- ? > 0 ---
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        J(vertices_i) = (1 - result_all{i}.p_gamma(2)) * (result_all{i}.gamma(2) > 0);
    end
    subplot(3,4,1 + (var_idx-1)*4);
    guide.Visualization.esi_plot(gca,J);
    colormap(pos_cmaps{var_idx});
    title([var_names{var_idx} ' \gamma > 0']);

    % --- ? < 0 ---
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        J(vertices_i) = (1 - result_all{i}.p_gamma(2)) * (result_all{i}.gamma(2) < 0);
    end
    %subplot(3,4,2 + (var_idx-1)*4);
    ax1 = subplot(3, 3, 1);
    guide.Visualization.esi_plot(ax1, J, colorLimits);
    %esi_plot_single;
    colormap(neg_cmaps{var_idx});
    title([var_names{var_idx} ' \gamma < 0']);

    % --- ß > 0 ---
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        J(vertices_i) = (1 - result_all{i}.p_beta(3)) * (result_all{i}.beta(3) > 0);
    end
    subplot(3,4,3 + (var_idx-1)*4);
    guide.Visualization.esi_plot(gca,J);
    colormap(pos_cmaps{var_idx});
    title([var_names{var_idx} ' \beta > 0']);

    % --- ß < 0 ---
    for i = 1:N_roi
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        J(vertices_i) = (1 - result_all{i}.p_beta(3)) * (result_all{i}.beta(3) < 0);
    end
    subplot(3,4,4 + (var_idx-1)*4);
    guide.Visualization.esi_plot(gca,J);
    colormap(neg_cmaps{var_idx});
    title([var_names{var_idx} ' \beta < 0']);
end

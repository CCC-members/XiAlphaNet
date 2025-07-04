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
num_groups = 6;       % Number of age groups (for age stratification)
mode = 1;             % 0 = amplitude plots, 1 = zero-inflation probability plots
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

ages = [];
All_Data = {}; 
index = 1;
for i=1:length(dataset.Participants)
    i
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

%% Cross Validation
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
        Xi_j_bin = Xi_j > threshold_Xi(j);
        PAF_j_bin = (PAF_j .* (Alpha_j > threshold_Alpha(j))) > threshold_PAF;
    else % For amplitude Distribution
        Alpha_j =  Alpha_j.*(Alpha_j > threshold_Alpha(j));
        Xi_j_bin = Xi_j.*(Xi_j > 0);
        PAF_j_bin = PAF_j.*((PAF_j .* (Alpha_j > threshold_Alpha(j))) > threshold_PAF);
    end

    AlphaAmp_all(:,j) = Alpha_j;
    XiAmp_all(:,j) = Xi_j_bin;
    PAF_all(:,j) = PAF_j_bin;
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age

if cross_index == 1

    % Define a grid of age values for estimation
    age_grid = linspace(min(ages), max(ages), 30);

    % Define the Epanechnikov kernel function
    kernel = @(u) (3/4)*(1 - u.^2) .* (abs(u) <= 1);

    % Define a range of bandwidths to search
    bandwidths = linspace(0.1, 100, 30);  % Adjust range and number of points as needed

    % Number of folds for cross-validation
    K = 10;

    % Initialize storage for cross-validation errors
    cv_error_PAF = zeros(length(bandwidths),1);
    cv_error_AlphaAmp = zeros(length(bandwidths),1);
    cv_error_XiAmp = zeros(length(bandwidths),1);

    % Prepare data for cross-validation
    data = struct();
    data.ages = ages;
    data.PAF = PAF_all(:);
    data.AlphaAmp = AlphaAmp_all(:);
    data.XiAmp = XiAmp_all(:);
    numSubjects = length(All_Data(1,:));
    % Create K-fold indices
    cv = cvpartition(numSubjects, 'KFold', K);

    % Iterate over each candidate bandwidth
    parfor b = 1:length(bandwidths)
        bandwidth = bandwidths(b);

        % Initialize temporary error accumulators
        temp_error_PAF = 0;
        temp_error_AlphaAmp = 0;
        temp_error_XiAmp = 0;

        for fold = 1:K
            % Define training and validation indices
            trainIdx = training(cv, fold);
            valIdx = test(cv, fold);

            % Extract training data
            ages_train = ages(trainIdx);
            PAF_train = PAF_all(:,trainIdx);
            AlphaAmp_train = AlphaAmp_all(:,trainIdx);
            XiAmp_train = XiAmp_all(:,trainIdx);

            % Extract validation data
            ages_val = ages(valIdx);
            PAF_val = PAF_all(:,valIdx);
            AlphaAmp_val = AlphaAmp_all(:,valIdx);
            XiAmp_val = XiAmp_all(:,valIdx);

            % Initialize predictions
            PAF_pred = zeros(size(PAF_val));
            AlphaAmp_pred = zeros(size(AlphaAmp_val));
            XiAmp_pred = zeros(size(XiAmp_val));

            % Perform kernel smoothing for each validation point
            for v = 1:length(ages_val)
                a = ages_val(v);
                u = (a - ages_train) / bandwidth;
                weights = kernel(u);

                % Normalize weights
                weights_sum = sum(weights);
                if weights_sum > 0
                    weights = weights / weights_sum;

                    % Predict for each measure
                    PAF_pred(:,v) = sum(PAF_train .* weights, 2);
                    AlphaAmp_pred(:,v) = sum(AlphaAmp_train .* weights, 2);
                    XiAmp_pred(:,v) = sum(XiAmp_train .* weights, 2);
                else
                    % If no weights, assign NaN or some default value
                    PAF_pred(:,v) = NaN;
                    AlphaAmp_pred(:,v) = NaN;
                    XiAmp_pred(:,v) = NaN;
                end
            end

            % Compute Mean Squared Error, ignoring NaNs
            valid = ~isnan(PAF_pred);
            temp_error_PAF = temp_error_PAF + mean((PAF_val(valid) - PAF_pred(valid)).^2);

            valid = ~isnan(AlphaAmp_pred);
            temp_error_AlphaAmp = temp_error_AlphaAmp + mean((AlphaAmp_val(valid) - AlphaAmp_pred(valid)).^2);

            valid = ~isnan(XiAmp_pred);
            temp_error_XiAmp = temp_error_XiAmp + mean((XiAmp_val(valid) - XiAmp_pred(valid)).^2);
        end

        % Average error across folds
        cv_error_PAF(b) = temp_error_PAF / K;
        cv_error_AlphaAmp(b) = temp_error_AlphaAmp / K;
        cv_error_XiAmp(b) = temp_error_XiAmp / K;

        fprintf('Bandwidth %.2f: CV Error PAF=%.4f, AlphaAmp=%.4f, XiAmp=%.4f\n', ...
            bandwidth, cv_error_PAF(b), cv_error_AlphaAmp(b), cv_error_XiAmp(b));
    end

    % Select the bandwidth with the minimum cross-validation error for each measure
    [~, idx_opt_PAF] = min(cv_error_PAF);
    [~, idx_opt_AlphaAmp] = min(cv_error_AlphaAmp);
    [~, idx_opt_XiAmp] = min(cv_error_XiAmp);

    h_t_PAF = bandwidths(idx_opt_PAF);
    h_t_AlphaAmp = bandwidths(idx_opt_AlphaAmp);
    h_t_XiAmp = bandwidths(idx_opt_XiAmp);

    fprintf('Optimal Bandwidths:\n');
    fprintf('PAF: %.2f\n', h_t_PAF);
    fprintf('AlphaAmp: %.2f\n', h_t_AlphaAmp);
    fprintf('XiAmp: %.2f\n', h_t_XiAmp);
else
    % Default Values 
    h_t_PAF = 100;
    h_t_AlphaAmp = 100;
    h_t_XiAmp = 100;
end
%% Kernel Density

% Define a grid of age values for estimation
age_grid = linspace(min(ages), max(ages), 40);  % 30 points for smoothness

% Define kernel parameters
% Assuming 'L' is a precomputed distance matrix of size [num_voxels x num_voxels]
L = parameters.Model.L;  % Distance matrix between voxels

% Define kernel functions
kernel_epanechnikov = @(u) (3/4)*(1 - u.^2) .* (abs(u) < 1);  % Epanechnikov Kernel

% Define bandwidths for spatial and temporal kernels
% Optimal Bandwidths (you can adjust these as needed)
h0 = mean(L(:));
h_s_PAF = h0;          % Spatial bandwidth for PAF 
h_s_AlphaAmp = h0;     % Spatial bandwidth for Alpha Amplitude 
h_s_XiAmp = h0;        % Spatial bandwidth for Xi Amplitude 


% Precompute spatial kernel matrices for each measure
% These are [num_voxels x num_voxels] matrices
spatial_kernel_PAF = kernel_epanechnikov(L / h_s_PAF);
spatial_kernel_AlphaAmp = kernel_epanechnikov(L / h_s_AlphaAmp);
spatial_kernel_XiAmp = kernel_epanechnikov(L / h_s_XiAmp);

% Initialize storage for the kernel estimates
num_voxels = size(PAF_all, 1);          % Number of voxels
num_age_points = length(age_grid);      % Number of age grid points

PAF_kernel = zeros(num_voxels, num_age_points);
AlphaAmp_kernel = zeros(num_voxels, num_age_points);
XiAmp_kernel = zeros(num_voxels, num_age_points);

% Replace zeros with NaN initially to handle cases with no data
PAF_kernel(:) = NaN;
AlphaAmp_kernel(:) = NaN;
XiAmp_kernel(:) = NaN;

% ----------------------------
% Spatio-Temporal Nadaraya-Watson Estimation
% ----------------------------

% Loop over each age in the age grid
for i = 1:num_age_points
    current_age = age_grid(i);
    
    % Compute temporal kernel weights for all subjects
    % [1 x num_subjects]
    u_t_PAF = (current_age - ages) / h_t_PAF;
    weights_t_PAF = kernel_epanechnikov(u_t_PAF);
    
    u_t_AlphaAmp = (current_age - ages) / h_t_AlphaAmp;
    weights_t_AlphaAmp = kernel_epanechnikov(u_t_AlphaAmp);
    
    u_t_XiAmp = (current_age - ages) / h_t_XiAmp;
    weights_t_XiAmp = kernel_epanechnikov(u_t_XiAmp);
    
    % Compute the sum of temporal weights for normalization
    sum_weights_t_PAF = sum(weights_t_PAF);
    sum_weights_t_AlphaAmp = sum(weights_t_AlphaAmp);
    sum_weights_t_XiAmp = sum(weights_t_XiAmp);
    
    % Handle cases where the sum of weights is zero to avoid division by zero
    if sum_weights_t_PAF == 0
        sum_weights_t_PAF = eps;  % A very small number
    end
    if sum_weights_t_AlphaAmp == 0
        sum_weights_t_AlphaAmp = eps;
    end
    if sum_weights_t_XiAmp == 0
        sum_weights_t_XiAmp = eps;
    end
    
    % Compute the weighted sums for each measure
    % [num_voxels x 1]
    weighted_PAF = PAF_all * weights_t_PAF';           % [num_voxels x 1]
    weighted_AlphaAmp = AlphaAmp_all * weights_t_AlphaAmp';
    weighted_XiAmp = XiAmp_all * weights_t_XiAmp';
    
    % Compute the Nadaraya-Watson estimates by applying spatial kernels
    % [num_voxels x 1] = [num_voxels x num_voxels] * [num_voxels x 1]
    sum_weighted_PAF = spatial_kernel_PAF * weighted_PAF;
    sum_weighted_AlphaAmp = spatial_kernel_AlphaAmp * weighted_AlphaAmp;
    sum_weighted_XiAmp = spatial_kernel_XiAmp * weighted_XiAmp;
    
    % Compute the normalization factors
    % [num_voxels x 1] = [num_voxels x num_voxels] * [1 x 1]
    normalization_PAF = sum(spatial_kernel_PAF,2) * sum(weights_t_PAF);
    normalization_AlphaAmp = sum(spatial_kernel_AlphaAmp,2) * sum(weights_t_AlphaAmp);
    normalization_XiAmp = sum(spatial_kernel_XiAmp,2) * sum(weights_t_XiAmp);
    
    % Calculate the density estimates
    PAF_estimate = sum_weighted_PAF ./ normalization_PAF;
    AlphaAmp_estimate = sum_weighted_AlphaAmp ./ normalization_AlphaAmp;
    XiAmp_estimate = sum_weighted_XiAmp ./ normalization_XiAmp;
    
    % Assign the estimates to the kernel storage matrices
    PAF_kernel(:,i) = PAF_estimate;
    AlphaAmp_kernel(:,i) = AlphaAmp_estimate;
    XiAmp_kernel(:,i) = XiAmp_estimate;
end
%% Marginalization 
age_group_edges = linspace(0,100,num_groups+1);

% Extract and preprocess ages
ages = cell2mat(All_Data(2,:));
ages(isnan(ages)) = mean(ages(~isnan(ages)));

% Assign each age to an age group
[~, ~, age_groups] = histcounts(ages, age_group_edges);

% Define the Epanechnikov kernel function
kernel = @(u) (3/4)*(1 - u.^2) .* (abs(u) < 1);


% Initialize storage for marginalized kernel estimates
PAF_marginalized = zeros(size(PAF_all, 1), num_groups);
AlphaAmp_marginalized = zeros(size(AlphaAmp_all, 1), num_groups);
XiAmp_marginalized = zeros(size(XiAmp_all, 1), num_groups);

% Define group centers
group_centers = (age_group_edges(1:end-1) + age_group_edges(2:end)) / 2;

% Loop through each entity and age group
parfor j = 1:size(PAF_all, 1)
    for g = 1:num_groups
        a = group_centers(g);
        group_idx = (age_groups == g);
        
        if any(group_idx)
            % PAF
            u_PAF = (a - ages(group_idx)) / h_t_PAF;
            weights_PAF = kernel(u_PAF);
            Z_PAF = sum(weights_PAF);
            if Z_PAF > 0
                PAF_marginalized(j, g) = sum(weights_PAF .* PAF_all(j, group_idx)) / Z_PAF;
            else
                PAF_marginalized(j, g) = NaN;
            end
            
            % AlphaAmp
            u_AlphaAmp = (a - ages(group_idx)) / h_t_AlphaAmp;
            weights_AlphaAmp = kernel(u_AlphaAmp);
            Z_AlphaAmp = sum(weights_AlphaAmp);
            if Z_AlphaAmp > 0
                AlphaAmp_marginalized(j, g) = sum(weights_AlphaAmp .* AlphaAmp_all(j, group_idx)) / Z_AlphaAmp;
            else
                AlphaAmp_marginalized(j, g) = NaN;
            end
            
            % XiAmp
            u_XiAmp = (a - ages(group_idx)) / h_t_XiAmp;
            weights_XiAmp = kernel(u_XiAmp);
            Z_XiAmp = sum(weights_XiAmp);
            if Z_XiAmp > 0
                XiAmp_marginalized(j, g) = sum(weights_XiAmp .* XiAmp_all(j, group_idx)) / Z_XiAmp;
            else
                XiAmp_marginalized(j, g) = NaN;
            end
        else
            PAF_marginalized(j, g) = NaN;
            AlphaAmp_marginalized(j, g) = NaN;
            XiAmp_marginalized(j, g) = NaN;
        end
    end
end

%% Plots
% === Step 1: Compute Average Estimates per Age Group ===
PAF_avg_intervals = PAF_marginalized;          % 1 x num_groups
AlphaAmp_avg_intervals = AlphaAmp_marginalized;% 1 x num_groups
XiAmp_avg_intervals = XiAmp_marginalized;      % 1 x num_groups
import guide.Visualization.esi_plot_single
% === Step 2: Define Age Intervals ===
age_intervals =  linspace(0,100,num_groups+1);% Adjust as needed

% === Step 3: Create the Plots ===
%figure('Position', [100, 100, 1500, 900]);  % Adjust figure size as needed

% === Plot PAF ===
for i = 1:max(num_groups-1,1)
    % subplot(3, num_groups, i);  % First row for PAF
    J_age_interval = PAF_avg_intervals(:,i);

    % Assign the current axis and set colormap
    %ax = gca;
    %colormap(ax, colormap_PAF);

    % Plot using custom function
    J = J_age_interval;
    p95 = prctile(J, 90);
    J(J > p95) = p95;
    esi_plot_single;  % Pass J as an argument
   % colormap('parula')
    % Set title and labels
    title(sprintf('PAF vs Age %.0f - %.0f', age_intervals(i), age_intervals(i+1)));
    ylabel('hot');
    %ylim([0, max_PAF * 1.1]);
    % xlabel('');  % Remove x-label for clarity
end

% === Plot AlphaAmp ===
for i = 1:max(num_groups-1,1)
    %subplot(3, num_groups, i + num_groups);  % Second row for AlphaAmp
    J_age_interval = AlphaAmp_avg_intervals(:,i);

    % Assign the current axis and set colormap
    %ax = gca;
    %colormap(ax, colormap_AlphaAmp);

    % Plot using custom function
    J = J_age_interval;
    p95 = prctile(J, 90);
    J(J > p95) = p95;
    esi_plot_single;  % Pass J as an argument
    %colormap('hot')
     %   colormap('parula')

    % Set title and labels
    title(sprintf('Alpha Amp vs Age %.0f - %.0f', age_intervals(i), age_intervals(i+1)));
    ylabel('Alpha Amp');
    %     ylim([0, max_AlphaAmp * 1.1]);
    %     xlabel('');  % Remove x-label for clarity
end

% === Plot XiAmp ===
for i = 1:max(num_groups-1,1)
    %subplot(3, num_groups, i + 2*num_groups);  % Third row for XiAmp
    J_age_interval = XiAmp_avg_intervals(:,i);

    % Assign the current axis and set colormap
    %     ax = gca;
    %     colormap(ax, colormap_XiAmp);

    % Plot using custom function
    J=J_age_interval;
    % Cap values above the 95th percentile
    p95 = prctile(J, 90);
    J(J > p95) = p95;
   esi_plot_single;  % Pass J as an argument
    %colormap('hot')   
   % colormap('parula')

    % Set title and labels
    title(sprintf('Xi Amp vs Age %.0f - %.0f', age_intervals(i), age_intervals(i+1)));
    ylabel('Xi Amp');

end




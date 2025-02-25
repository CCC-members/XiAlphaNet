clc;
clear all;
% 
% % Paths
% DataPath = 'Data\Scalp_Density_Matrix';
% modelParametersPath = 'Data\Model_Parameters';
% 
% % Subfolders within the main folder
% subFolders = {'Control2'};
% 
% % Load parameters
% load('Data\Model_Parameters\parameters.mat');
% Ne = parameters.Dimensions.Ne;
% Nr = parameters.Dimensions.Nr;
% Nv = parameters.Dimensions.Nv;
% Nw = parameters.Dimensions.Nw;
% %
% All_Data = {};  % Initialize a cell array to store data
% x_avg=zeros(7*Nv+1,1899);
% % Loop through the subfolders and process the data
% for k = 1:length(subFolders)
%     folderPath_X = fullfile(modelParametersPath, subFolders{k});
%     matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
% 
%     for j = 1:length(matFiles_X)
%         filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
%         data_X = load(filePath_X);
% 
%         % Extract and store the data
%         %x_avg(:,j) =  data_X.x.Solution;
%         All_Data{1,j} = data_X.x.Solution;
%         All_Data{2,j} = data_X.x.Age;
%     end
% end

age_min = 0;%age_range(1);
age_max = 100;%age_range(2);
dataset = jsondecode(fileread('/home/ronaldo/Documents/dev/Data/Results/XIALPHANET.json'));
parameters = load('/home/ronaldo/Documents/dev/Data/Results/structural/parameters.mat');
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
       % if (isequal(Process,'Alpha'))
       %     spec_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
       % elseif (isequal(Process,'Xi'))
       %     spec_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
       % end
       % if isequal(variable,'Power')
       %      All_Data{1,i} = spec_process.Power;
       % elseif (isequal(variable,'Width'))
       %      All_Data{1,i} = spec_process.Width;
       % elseif (isequal(variable,'Exponent'))
       %      All_Data{1,i} = spec_process.Exponent;
       % elseif (isequal(variable,'PAF')) && (isequal(Process,'Alpha'))
       %      All_Data{1,i} = spec_process.PAF;
       % end
    end
end


%% Simple Counting 
% Initialize storage for Peak Alpha Frequency (a(:,4)), Amplitude of the Alpha (a(:,1)), and Amplitude of Xi
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);                                % Store Amplitude of the Xi
end


% Define threshold for each measure use set_threshold_em for recalculate this value
threshold_PAF = 7;             % Example threshold for PAF
% threshold_AlphaAmp = 0.19;   % Example threshold for Alpha Amplitude
% threshold_XiAmp = 0.19;      % Example threshold for Xi Amplitude

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold_AlphaAmp   =  prctile(AlphaAmp_all(:,j),75);
    PAF_all(:,j) = PAF_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = XiAmp_all(:,j);         % Store Amplitude of the Xi
end


% Apply threshold and create binary masks for each measure
for j = 1:length(All_Data(1,:))
    fprintf('Processing subject %d\n', j);
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold_AlphaAmp   =  prctile(AlphaAmp_all(:,j),75);
    PAF_all(:,j) = PAF_all(:,j) > threshold_PAF;             % Apply threshold for PAF
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j) > threshold_AlphaAmp;   % Apply threshold for Alpha Amplitude
    XiAmp_all(:,j) =  XiAmp_all(:,j) > threshold_AlphaAmp;         % Apply threshold for Xi Amplitude
end


% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages

% Define age intervals
age_intervals = linspace(age_min, age_max, floor((age_max-age_min)/20)+1);  % 5 intervals, 6 edges

% Initialize storage for the averaged values
PAF_avg_intervals = zeros(size(PAF_all, 1), length(age_intervals)-1);
AlphaAmp_avg_intervals = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
XiAmp_avg_intervals = zeros(size(XiAmp_all, 1), length(age_intervals)-1);

% Average values for each age interval
for i = 1:length(age_intervals)-1
    age_mask = (ages >= age_intervals(i)) & (ages < age_intervals(i+1));
    
    % For PAF
    f_va_PAF = sum(PAF_all(:, age_mask), 2);
    nf_va_PAF = sum(age_mask);  % Number of subjects in this age group
    PAF_avg_intervals(:, i) = f_va_PAF / (1+nf_va_PAF);
    
    % For Alpha Amplitude
    f_va_AlphaAmp = sum(AlphaAmp_all(:, age_mask), 2);
    nf_va_AlphaAmp = sum(age_mask);  % Number of subjects in this age group
    AlphaAmp_avg_intervals(:, i) = f_va_AlphaAmp / (1+nf_va_AlphaAmp);
    
    % For Xi Amplitude
    f_va_XiAmp = sum(XiAmp_all(:, age_mask), 2);
    nf_va_XiAmp = sum(age_mask);  % Number of subjects in this age group
    XiAmp_avg_intervals(:, i) = f_va_XiAmp / (1+nf_va_XiAmp);
end
%%
% Create a single large figure for all subplots
%figure('Position', [100, 100, 1500, 900]);  % Adjust figure size for better visibility

% Plotting for Peak Alpha Frequency (PAF) in the first row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*PAF_avg_intervals(:,i) ;
    J = J_age_interval;
   % subplot(3, 5, i);  % 3 rows, 5 columns, first row for PAF
    %esi_plot(gca, J_age_interval, [0, max(PAF_avg_intervals(:))]);
    esi_plot_single
    title(sprintf('Alpha Peak Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Plotting for Alpha Amplitude in the second row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*AlphaAmp_avg_intervals(:,i) ;
    
    %subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
    %esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
    J = J_age_interval;
    esi_plot_single
    % colormap("magma");
    title(sprintf('Alpha Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Plotting for Xi Amplitude in the third row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*XiAmp_avg_intervals(:,i) ;
    
    %subplot(3, 5, i + 10);  % 3 rows, 5 columns, third row for Xi Amplitude
    %esi_plot(gca, J_age_interval, [0, max(XiAmp_avg_intervals(:))]);
    %colormap("hot");
    J = J_age_interval;
    esi_plot_single
    title(sprintf('Xi Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Adjust layout and appearance
sgtitle('Lifespan Mapping of Probability Distribution');  % Overall title

% Create a single large figure for all subplots
figure('Position', [100, 100, 1500, 900]);  % Adjust figure size for better visibility

% % Plotting for Peak Alpha Frequency (PAF) in the first row
% for i = 1:length(age_intervals)-1
%     J_age_interval = 1.*(PAF_avg_intervals(:,i) > 0.90);
%     J = J_age_interval;
%    % subplot(3, 5, i);  % 3 rows, 5 columns, first row for PAF
%    % esi_plot(gca, J_age_interval, [0, max(PAF_avg_intervals(:))]);
%    esi_plot_single
%     %colormap("hot");
%     title(sprintf('Alpha Peak Freq vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
% end

% Plotting for Alpha Amplitude in the second row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*AlphaAmp_avg_intervals(:,i) > 0.99;
    J = J_age_interval;
    %subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
    %esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
    esi_plot_single
    title(sprintf('Alpha Peak Freq vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
    %colormap("hot");
end

% Plotting for Xi Amplitude in the third row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*(XiAmp_avg_intervals(:,i) > 0.95);
    J = J_age_interval;
    %subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
    %esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
    esi_plot_single
    title(sprintf('Xi Amplitude vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Adjust layout and appearance
sgtitle('Lifespan Mapping of High-Probability Voxel Activity (p>0.5)');  % Overall title

% Create a single large figure for all subplots
figure('Position', [100, 100, 1500, 900]);  % Adjust figure size for better visibility


%% Kernel Density
% Initialize storage for Peak Alpha Frequency (PAF), Amplitude of the Alpha, and Amplitude of Xi
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);                                % Store Amplitude of the Xi
end

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold_AlphaAmp   =  prctile(AlphaAmp_all(:,j),25);
    PAF_all(:,j) = PAF_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = XiAmp_all(:,j);         % Store Amplitude of the Xi
end


% Apply threshold and create binary masks for each measure
for j = 1:length(All_Data(1,:))
    fprintf('Processing subject %d\n', j);
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold_AlphaAmp   =  prctile(AlphaAmp_all(:,j),25);
    PAF_all(:,j) = PAF_all(:,j) > threshold_PAF;             % Apply threshold for PAF
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j) > threshold_AlphaAmp;   % Apply threshold for Alpha Amplitude
    XiAmp_all(:,j) =  XiAmp_all(:,j) > threshold_AlphaAmp;         % Apply threshold for Xi Amplitude
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age

% Define a grid of age values for estimation
age_grid = linspace(min(ages), max(ages), floor((max(ages)-min(ages))/100*30)+1);  % 30 points for smoothness

% Define kernel parameters
% Assuming 'L' is a precomputed distance matrix of size [num_voxels x num_voxels]
L = parameters.Model.L;  % Distance matrix between voxels

% Define kernel functions
kernel_epanechnikov = @(u) (3/4)*(1 - u.^2) .* (abs(u) < 1);  % Epanechnikov Kernel

% Define bandwidths for spatial and temporal kernels
% Optimal Bandwidths (you can adjust these as needed)
h_s_PAF = 1;          % Spatial bandwidth for PAF 
h_s_AlphaAmp = 100;     % Spatial bandwidth for Alpha Amplitude 
h_s_XiAmp = 100;        % Spatial bandwidth for Xi Amplitude 

% Temporal bandwidths (assuming they are the same as spatial for simplicity)
% You can define separate temporal bandwidths if needed
% Alternatively, define h_t differently for each measure
h_t_PAF = 1.73;
h_t_AlphaAmp = 1.73;
h_t_XiAmp = 1.73;

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
age_group_edges = age_min:(floor((age_max-age_min)/20)+1):age_max;
num_groups = length(age_group_edges) - 1;

% Extract and preprocess ages
ages = cell2mat(All_Data(2,:));
ages(isnan(ages)) = mean(ages(~isnan(ages)));

% Assign each age to an age group
[~, ~, age_groups] = histcounts(ages, age_group_edges);

% Define the Epanechnikov kernel function
kernel = @(u) (3/4)*(1 - u.^2) .* (abs(u) < 1);

% Define optimal bandwidths
h_t_PAF = 1.73;
h_t_AlphaAmp = 1.73;
h_t_XiAmp = 1.73;

% Initialize storage for marginalized kernel estimates
PAF_marginalized = zeros(size(PAF_all, 1), num_groups);
AlphaAmp_marginalized = zeros(size(AlphaAmp_all, 1), num_groups);
XiAmp_marginalized = zeros(size(XiAmp_all, 1), num_groups);

% Define group centers
group_centers = (age_group_edges(1:end-1) + age_group_edges(2:end)) / 2;

% Loop through each entity and age group
for j = 1:size(PAF_all, 1)
    for g = 1:num_groups
        a = group_centers(g);
        group_idx = (age_groups == g);
        
        if any(group_idx)
            % PAF
            u_PAF = (a - ages(group_idx)) / h_s_PAF;
            weights_PAF = kernel(u_PAF);
            Z_PAF = sum(weights_PAF);
            if Z_PAF > 0
                PAF_marginalized(j, g) = sum(weights_PAF .* PAF_all(j, group_idx)) / Z_PAF;
            else
                PAF_marginalized(j, g) = NaN;
            end
            
            % AlphaAmp
            u_AlphaAmp = (a - ages(group_idx)) / h_s_AlphaAmp;
            weights_AlphaAmp = kernel(u_AlphaAmp);
            Z_AlphaAmp = sum(weights_AlphaAmp);
            if Z_AlphaAmp > 0
                AlphaAmp_marginalized(j, g) = sum(weights_AlphaAmp .* AlphaAmp_all(j, group_idx)) / Z_AlphaAmp;
            else
                AlphaAmp_marginalized(j, g) = NaN;
            end
            
            % XiAmp
            u_XiAmp = (a - ages(group_idx)) / h_s_XiAmp;
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

% === Step 2: Define Age Intervals ===
age_intervals =  age_min:(floor((age_max-age_min)/20)+1):age_max;% Adjust as needed
num_groups = length(age_intervals) - 1;

% === Step 3: Create the Plots ===
%figure('Position', [100, 100, 1500, 900]);  % Adjust figure size as needed

%%
% === Plot PAF ===
for i = 1:num_groups-1
   % subplot(3, num_groups, i);  % First row for PAF
    J_age_interval = PAF_avg_intervals(:,i);
    
    % Assign the current axis and set colormap
    %ax = gca;
    %colormap(ax, colormap_PAF);
    
    % Plot using custom function
    J = J_age_interval;
    esi_plot_single;  % Pass J as an argument
    
    % Set title and labels
    title(sprintf('PAF vs Age %.0f - %.0f', age_intervals(i), age_intervals(i+1)));
    ylabel('PAF');
    %ylim([0, max_PAF * 1.1]);
   % xlabel('');  % Remove x-label for clarity
end

% === Plot AlphaAmp ===
for i = 1:num_groups
    %subplot(3, num_groups, i + num_groups);  % Second row for AlphaAmp
    J_age_interval = AlphaAmp_avg_intervals(:,i);
    
    % Assign the current axis and set colormap
    %ax = gca;
    %colormap(ax, colormap_AlphaAmp);
    
    % Plot using custom function
    J = J_age_interval;
    esi_plot_single;  % Pass J as an argument
    
    % Set title and labels
    title(sprintf('Alpha Amp vs Age %.0f - %.0f', age_intervals(i), age_intervals(i+1)));
    ylabel('Alpha Amp');
%     ylim([0, max_AlphaAmp * 1.1]);
%     xlabel('');  % Remove x-label for clarity
end

% === Plot XiAmp ===
for i = 1:num_groups
    %subplot(3, num_groups, i + 2*num_groups);  % Third row for XiAmp
    J_age_interval = XiAmp_avg_intervals(:,i);
    
    % Assign the current axis and set colormap
%     ax = gca;
%     colormap(ax, colormap_XiAmp);
    
    % Plot using custom function
    J=J_age_interval;
    esi_plot_single;  % Pass J as an argument
    
    % Set title and labels
    title(sprintf('Xi Amp vs Age %.0f - %.0f', age_intervals(i), age_intervals(i+1)));
    ylabel('Xi Amp');
    
end

% === Add Overall Title ===
sgtitle('Lifespan Mapping of Probability Distribution');  % Overall title

% === Optional: Adjust Layout ===
% MATLAB does not have a built-in 'tight_layout' like Python's matplotlib.
% To improve spacing, you can adjust subplot parameters manually.
% Alternatively, use the 'sgtitle' to minimize overlap.

% Example: Adjust subplot positions
% for idx = 1:(3*num_groups)
%     subplot_handles(idx) = subplot(3, num_groups, idx);
%     set(subplot_handles(idx), 'LooseInset', get(gca, 'TightInset'));
% end



%% Cross Validation
% Initialize storage for Peak Alpha Frequency (PAF), Amplitude of the Alpha, and Amplitude of Xi
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);                                % Store Amplitude of the Xi
end

% Calculate sparse index (optional, can be used for any of the measures)
sparse_index_PAF = sum(PAF_all(:) == 0) / numel(PAF_all) * 100;
sparse_index_AlphaAmp = sum(AlphaAmp_all(:) == 0) / numel(AlphaAmp_all) * 100;
sparse_index_XiAmp = sum(XiAmp_all(:) == 0) / numel(XiAmp_all) * 100;

% Define threshold for each measure use set_threshold_em for recalculate this value
threshold_PAF = 7;             % Example threshold for PAF
threshold_AlphaAmp = 0.19;   % Example threshold for Alpha Amplitude
threshold_XiAmp = 0.19;      % Example threshold for Xi Amplitude

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = PAF_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = XiAmp_all(:,j);         % Store Amplitude of the Xi
end


% Apply threshold and create binary masks for each measure
for j = 1:length(All_Data(1,:))
    fprintf('Processing subject %d\n', j);
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = PAF_all(:,j) > threshold_PAF;             % Apply threshold for PAF
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j) > threshold_AlphaAmp;   % Apply threshold for Alpha Amplitude
    XiAmp_all(:,j) =  XiAmp_all(:,j) > threshold_XiAmp;         % Apply threshold for Xi Amplitude
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age

% Define a grid of age values for estimation
age_grid = linspace(min(ages), max(ages), 40);  % 40 points for grid

% Define the Epanechnikov kernel function
kernel = @(u) (3/4)*(1 - u.^2) .* (abs(u) <= 1);

% Define a range of bandwidths to search
bandwidths = linspace(0.1, 5, 10);  % Adjust range and number of points as needed

% Number of folds for cross-validation
K = 5;

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
numSubjects = 1965;
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

h_s_PAF = bandwidths(idx_opt_PAF);
h_s_AlphaAmp = bandwidths(idx_opt_AlphaAmp);
h_s_XiAmp = bandwidths(idx_opt_XiAmp);

fprintf('Optimal Bandwidths:\n');
fprintf('PAF: %.2f\n', h_s_PAF);
fprintf('AlphaAmp: %.2f\n', h_s_AlphaAmp);
fprintf('XiAmp: %.2f\n', h_s_XiAmp);

%Optimal Bandwidths:
%PAF: 4.46
%AlphaAmp: 3.37
%XiAmp: 1.73

%%
% ----------------------------
% Spatio-Temporal Density Estimation using Nadaraya-Watson Estimator
% ----------------------------

% Initialize storage for Peak Alpha Frequency (PAF), Amplitude of the Alpha, and Amplitude of Xi
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4) .* (a(:,1) > 1.8632);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1) .* (a(:,1) > 1.8632);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);                                % Store Amplitude of the Xi
end

% Calculate sparse index (optional, can be used for any of the measures)
sparse_index_PAF = sum(PAF_all(:) == 0) / numel(PAF_all) * 100;
sparse_index_AlphaAmp = sum(AlphaAmp_all(:) == 0) / numel(AlphaAmp_all) * 100;
sparse_index_XiAmp = sum(XiAmp_all(:) == 0) / numel(XiAmp_all) * 100;

% Define threshold for each measure use set_threshold_em for recalculate this value
threshold_PAF = 7;             % Example threshold for PAF
threshold_AlphaAmp = 1.8632;   % Example threshold for Alpha Amplitude
threshold_XiAmp = 1.8632;      % Example threshold for Xi Amplitude

% Apply threshold and create binary masks for each measure
for j = 1:length(All_Data(1,:))
    fprintf('Processing subject %d\n', j);
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4) > threshold_PAF;             % Apply threshold for PAF
    AlphaAmp_all(:,j) = a(:,1) > threshold_AlphaAmp;   % Apply threshold for Alpha Amplitude
    XiAmp_all(:,j) = e(:,1) > threshold_XiAmp;         % Apply threshold for Xi Amplitude
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age

% Define a grid of age values for estimation
age_grid = linspace(min(ages), max(ages), 30);  % 30 points for smoothness

% Define kernel parameters
% Assuming 'L' is a precomputed distance matrix of size [num_voxels x num_voxels]
L = parameters.Model.L;  % Distance matrix between voxels

% Define kernel functions
kernel_epanechnikov = @(u) (3/4)*(1 - u.^2) .* (abs(u) < 1);  % Epanechnikov Kernel

% Define bandwidths for spatial and temporal kernels
% Optimal Bandwidths (you can adjust these as needed)
h_s_PAF = 0.1;          % Spatial bandwidth for PAF 
h_s_AlphaAmp = 0.100;     % Spatial bandwidth for Alpha Amplitude 
h_s_XiAmp = 0.100;        % Spatial bandwidth for Xi Amplitude 

% Temporal bandwidths (assuming they are the same as spatial for simplicity)
% You can define separate temporal bandwidths if needed
% Alternatively, define h_t differently for each measure
h_t_PAF = 4.46;
h_t_AlphaAmp = 3.37;
h_t_XiAmp = 1.73;

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





function plot_peak_freq_vs_age_vs_region(region, hemisphere)
% Initialize the vertices array
vertices = [];

% Check if both region and hemisphere inputs are provided
if isempty(region) || isempty(hemisphere)
    error('Please specify both region and hemisphere (L or R).');
end

% Define region lists based on the HCP MMP1 Atlas
switch region

    case 'Prefrontal'
        % Prefrontal Cortex Regions
        if hemisphere == 'L'
            regions = {'L_10d_ROI L', 'L_10pp_ROI L', 'L_10r_ROI L', 'L_a10p_ROI L'};
        elseif hemisphere == 'R'
            regions = {'R_10d_ROI R', 'R_10pp_ROI R', 'R_10r_ROI R', 'R_a10p_ROI R'};
        end
   
    case 'Temporal'
        % Temporal Lobe Middle Areas
        if hemisphere == 'L'
            regions = {'L_TE1a_ROI L', 'L_TE1p_ROI L', 'L_TE2a_ROI L', 'L_TE2p_ROI L'};
        elseif hemisphere == 'R'
            regions = {'R_TE1a_ROI R', 'R_TE1p_ROI R', 'R_TE2a_ROI R','R_TE2p_ROI R'};
        end

    case 'Parietal'
        % Superior Parietal Lobe Regions
        if hemisphere == 'L'
            regions = {'L_7AL_ROI L', 'L_7PC_ROI L', 'L_7PL_ROI L', 'L_7Pm_ROI L'};
        elseif hemisphere == 'R'
            regions = {'R_7AL_ROI R', 'R_7PC_ROI R', 'R_7PL_ROI R', 'R_7Pm_ROI R'};
        end

    case 'Occipital'
        % Occipital Lobe Regions
        if hemisphere == 'L'
            regions = {'L_V1_ROI L', 'L_V2_ROI L', 'L_V3A_ROI L', 'L_V3B_ROI L'};
        elseif hemisphere == 'R'
            regions = {'R_V1_ROI R', 'R_V2_ROI R', 'R_V3A_ROI R', 'R_V3B_ROI R'};
        end

    otherwise
        error('Invalid region specified. Choose from: Prefrontal, Temporal, Parietal, Occipital.');
        
end

% Extract vertices for the specified regions
for i = 1:length(regions)
    vertices = [vertices, getRegionVertices(regions{i})];
end

% Define paths
%DataPath = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';
%dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control'};

% Load necessary parameters
load('Data\Model_Parameters\parameters.mat');

%load('Data\RinvT\RinvT_G.mat', 'RinvT');
%J=zeros(8003,1);
%J(vertices) = 1;
%esi_plot_single;
% Initialize cell array
All_Data = {};
%%
% Process data
for k = 1:length(subFolders)
    % Full path to the subfolder
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
    
    % Second pass to process the data with the new threshold
    index = 1; % Reset index
    for j = randperm(length(matFiles_X),length(matFiles_X))%1:length(matFiles_X)%
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        data = load(filePath_X);
    
        [e, a, s2] = x2v(data.x.Solution);
    
        % Apply the calculated threshold to find vertices
        threshold = prctile(a(:,1),5);
        % 2. Find all indices where amplitude > threshold
         %activeIndices = find(a(:,1) > 0);
        
        % 3. Find the intersection between 'vertices' and 'activeIndices'
        %vert2 = intersect(vertices, activeIndices);
        a(:,4)=a(:,4).*(a(:,1)>threshold);
        J = a(vertices, 4);%.*(a(vertices,1)>threshold);
        % Exclude possible NaNs
        J = J(~isnan(J));
    
        % Store in All_Data
        All_Data{1, index} = J; % Peak alpha frequencies
        All_Data{2, index} = data.x.Age; % Age of the subject
        index = index + 1;
    end

end

% Perform a robust measure of central tendency
% For example, use the median or a trimmed mean to avoid the influence of outliers
% Extract ages and calculate mean voxel activation per subject
% Initialize arrays to collect voxel activations, ages, and voxel positions for all subjects
% Collect voxel activations, ages, and voxel positions
%%
all_voxel_activations = [];
all_ages = [];
voxel_positions = [];

% Loop through each subject to collect data
for j = 1:length(All_Data(1, :))
    voxel_data = All_Data{1, j};  % Voxel activation data for subject j
    age = All_Data{2, j};         % Age for subject j
    
    % Get number of voxels and their IDs
    num_voxels = numel(voxel_data);
    voxel_ids = (1:num_voxels)';  % Unique IDs for each voxel
    
    % Append to overall data
    all_voxel_activations = [all_voxel_activations; voxel_data(:)];
    all_ages = [all_ages; repmat(age, num_voxels, 1)];
    voxel_positions = [voxel_positions; voxel_ids]; % Assuming voxel IDs are random effects
end

% Define the age range for evaluation points (e.g., from min to max age)
age_range = linspace(1, 90, 100)';
%%
% Call the function to fit the zero-inflated model with random effects
model_output = fit_zero_inflation_random_effects(all_ages, all_voxel_activations, voxel_positions, age_range);

%%
x=age_range;


cond_pred = model_output.predictions.cond;  % Conditional model predictions
cond_se = model_output.predictions.cond_se;  % Conditional model standard errors
zi_pred = model_output.predictions.zi;  % Zero-inflation model predictions (probabilities)
zi_se =  model_output.predictions.zi_se;  % Zero-inflation model standard errors

%Replace NaN values with the mean of the non-NaN values
cond_se(isnan(cond_se)) = nanmean(cond_se(~isnan(cond_se)));  % Replace NaNs with mean in cond_se
zi_se(isnan(zi_se)) = nanmean(zi_se(~isnan(zi_se)));  % Replace NaNs with mean in zi_se

if length(cond_se)~= cond_pred
    cond_se = min(std(all_voxel_activations-mean(cond_se)),0.1)*ones(1,100);
    zi_se  = 0.01*ones(1,100);
   % disp('Not Zero')
end
cond_se = movmean(cond_se,20);
zi_se = movmean(zi_se,20);

% Calculate the upper and lower bounds for the final prediction
% The final prediction is a mixture of the conditional model prediction and zero-inflation probability

final_pred = (1 - zi_pred) .* cond_pred;  % Final prediction: conditional model prediction weighted by (1 - zero-inflation probability)

% Upper and lower bounds for the final prediction
if mean(cond_se)< 10^(-2)
    sigma1 = 10;
else
    sigma1 = 2;
end

cond_upper = (cond_pred + sigma1*cond_se);  % Upper bound for final prediction
cond_lower =  (cond_pred - sigma1*cond_se);  % Lower bound for final prediction

final_upper = (1 - zi_pred) .* (cond_pred + sigma1*cond_se);  % Upper bound for final prediction
final_lower =  (1 - zi_pred) .* (cond_pred - sigma1*cond_se);  % Lower bound for final prediction

if mean(cond_se)< 10^(-2) 
    sigma2 = 10;
else
    sigma2 = 2;
end

% Upper and lower bounds for the final prediction
zi_final_upper =  zi_pred + sigma2*reshape(zi_se,size(zi_pred));  % Upper bound for final prediction
zi_final_lower =  zi_pred - sigma2*reshape(zi_se,size(zi_pred));  % Lower bound for final prediction

%cond_pred = final_pred;
%%
% ------------------------------
% Plotting Final Predictions with Error Bands
% ------------------------------

% ------------------------------
% Plotting Final Predictions with Error Bands
% Zero Inflation Model Results
% ------------------------------

% Define the main title for the entire figure
mainTitle = region;  % Editable main title

% Create a figure (optional, but recommended for clarity)
figure;

% ------------------------------
% Subplot 1: Conditional Prediction
% ------------------------------
if hemisphere == 'L'
subplot(3,1,1);  % Select the first subplot in a 2x1 grid
hold on;         % Hold the current axes for multiple plots

% Ensure that x, final_upper, and final_lower are column vectors
x = x(:);                     % Convert x to a column vector
cond_upper = cond_upper(:); % Convert final_upper to a column vector
cond_lower = cond_lower(:); % Convert final_lower to a column vector

% Plot the shaded error bands (filled area)
fill([x; flipud(x)], [cond_upper; flipud(cond_lower)], 'b', ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the final prediction (central curve)
plot(x, cond_pred, 'b', 'LineWidth', 2);

% Label the x and y axes with bold font
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('APF (Hz)', 'FontWeight', 'bold');

% Add title to the subplot with bold font
title('Conditional Prediction (E[Y|Y>0])', 'FontWeight', 'bold');

% Set axis limits based on data ranges with some padding
axis([min(x), max(x), min(cond_pred)-0.5, max(cond_pred)+0.5 ]);

% Set the background color of the axes to white
set(gca, 'Color', 'w');

% Enable grid for better readability
grid on;

% Make axis numbers bold
set(gca, 'FontWeight', 'bold');

% Release the hold on the current axes
hold off;

% ------------------------------
% Subplot 2: Zero-Inflation Prediction
% ------------------------------
subplot(3,1,2);  % Select the second subplot in a 2x1 grid
hold on;         % Hold the current axes for multiple plots

% Ensure that x, zi_final_upper, and zi_final_lower are column vectors
x = x(:);                          % Convert x to a column vector
zi_final_upper = zi_final_upper(:);% Convert zi_final_upper to a column vector
zi_final_lower = zi_final_lower(:);% Convert zi_final_lower to a column vector

% Plot the shaded error bands (filled area)
fill([x; flipud(x)], [zi_final_upper; flipud(zi_final_lower)], 'b', ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the final prediction (central curve)
plot(x, zi_pred, 'b', 'LineWidth', 2);
plot(x,0.5*ones(size(x)),'g--', 'LineWidth', 2)

% Label the x and y axes with bold font
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('Probability', 'FontWeight', 'bold');

% Add title to the subplot with bold font
title('Zero-Inflation Prediction  (P[Y=0])', 'FontWeight', 'bold');

% Set axis limits based on data ranges with appropriate bounds
axis([min(x), max(x), 0, 1]); 

% Set the background color of the axes to white
set(gca, 'Color', 'w');

% Enable grid for better readability
grid on;

% Make axis numbers bold
set(gca, 'FontWeight', 'bold');

% Release the hold on the current axes
hold off;
subplot(3,1,3);  % Select the first subplot in a 2x1 grid
hold on;         % Hold the current axes for multiple plots

% Ensure that x, final_upper, and final_lower are column vectors
x = x(:);                     % Convert x to a column vector
final_upper = final_upper(:); % Convert final_upper to a column vector
final_lower = final_lower(:); % Convert final_lower to a column vector

% Plot the shaded error bands (filled area)
fill([x; flipud(x)], [final_upper; flipud(final_lower)], 'b', ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the final prediction (central curve)
plot(x, final_pred, 'b', 'LineWidth', 2);

% Label the x and y axes with bold font
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('APF (Hz)', 'FontWeight', 'bold');

% Add title to the subplot with bold font
title('Conditional Prediction (E[Y|Y>0])', 'FontWeight', 'bold');

% Set axis limits based on data ranges with some padding
axis([min(x), max(x), min(final_pred)-0.5, max(final_pred)+0.5 ]);

% Set the background color of the axes to white
set(gca, 'Color', 'w');

% Enable grid for better readability
grid on;

% Make axis numbers bold
set(gca, 'FontWeight', 'bold');

% Release the hold on the current axes
hold off;

% ------------------------------
% Add Main Title to the Entire Figure
% ------------------------------
% Check if 'sgtitle' is available (introduced in R2018b)
if exist('sgtitle', 'file')
    sgtitle(mainTitle, 'FontWeight', 'bold', 'FontSize', 14);
else
    % For older MATLAB versions, use a workaround
    % Add an invisible subplot to hold the main title
    subplot(2,1,1);
    title('');
    subplot(2,1,2);
    title('');
    annotation('textbox', [0.5, 0.95, 0, 0], ...
               'String', mainTitle, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'FontWeight', 'bold', ...
               'FontSize', 14);
end
else
  subplot(3,1,1);  % Select the first subplot in a 2x1 grid
hold on;         % Hold the current axes for multiple plots

% Ensure that x, final_upper, and final_lower are column vectors
x = x(:);                     % Convert x to a column vector
cond_upper = cond_upper(:); % Convert final_upper to a column vector
cond_lower = cond_lower(:); % Convert final_lower to a column vector

% Plot the shaded error bands (filled area)
fill([x; flipud(x)], [cond_upper; flipud(cond_lower)], 'r', ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the final prediction (central curve)
plot(x, cond_pred, 'r', 'LineWidth', 2);

% Label the x and y axes with bold font
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('APF (Hz)', 'FontWeight', 'bold');

% Add title to the subplot with bold font
title('Prediction ((1-P[Y=0])\times E[Y|Y>0])', 'FontWeight', 'bold');

% Set axis limits based on data ranges with some padding
axis([min(x), max(x),  min(cond_pred)-0.5, max(cond_pred)+0.5 ]);

% Set the background color of the axes to white
set(gca, 'Color', 'w');

% Enable grid for better readability
grid on;

% Make axis numbers bold
set(gca, 'FontWeight', 'bold');

% Release the hold on the current axes
hold off;

% ------------------------------
% Subplot 2: Zero-Inflation Prediction
% ------------------------------
subplot(3,1,2);  % Select the second subplot in a 2x1 grid
hold on;         % Hold the current axes for multiple plots

% Ensure that x, zi_final_upper, and zi_final_lower are column vectors
x = x(:);                          % Convert x to a column vector
zi_final_upper = zi_final_upper(:);% Convert zi_final_upper to a column vector
zi_final_lower = zi_final_lower(:);% Convert zi_final_lower to a column vector

% Plot the shaded error bands (filled area)
fill([x; flipud(x)], [zi_final_upper; flipud(zi_final_lower)], 'r', ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the final prediction (central curve)
plot(x, zi_pred, 'r', 'LineWidth', 2);
plot(x,0.5*ones(size(x)),'g--', 'LineWidth', 2)
% Label the x and y axes with bold font
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('Probability', 'FontWeight', 'bold');

% Add title to the subplot with bold font
title('Zero-Inflation Prediction (P[Y=0])', 'FontWeight', 'bold');

% Set axis limits based on data ranges with appropriate bounds
axis([min(x), max(x), 0, 1]); 

% Set the background color of the axes to white
set(gca, 'Color', 'w');

% Enable grid for better readability
grid on;

% Make axis numbers bold
set(gca, 'FontWeight', 'bold');
subplot(3,1,3);  % Select the first subplot in a 2x1 grid
hold on;         % Hold the current axes for multiple plots

% Ensure that x, final_upper, and final_lower are column vectors
x = x(:);                     % Convert x to a column vector
final_upper = final_upper(:); % Convert final_upper to a column vector
final_lower = final_lower(:); % Convert final_lower to a column vector

% Plot the shaded error bands (filled area)
fill([x; flipud(x)], [final_upper; flipud(final_lower)], 'r', ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the final prediction (central curve)
plot(x, final_pred, 'r', 'LineWidth', 2);

% Label the x and y axes with bold font
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('APF (Hz)', 'FontWeight', 'bold');

% Add title to the subplot with bold font
title('Prediction ((1-P[Y=0])\times E[Y|Y>0])', 'FontWeight', 'bold');

% Set axis limits based on data ranges with some padding
axis([min(x), max(x), min(final_pred)-0.5, max(final_pred)+0.5 ]);

% Set the background color of the axes to white
set(gca, 'Color', 'w');

% Enable grid for better readability
grid on;

% Make axis numbers bold
set(gca, 'FontWeight', 'bold');

% Release the hold on the current axes
hold off;

% ------------------------------
% Add Main Title to the Entire Figure
% ------------------------------
% Check if 'sgtitle' is available (introduced in R2018b)
if exist('sgtitle', 'file')
    sgtitle(mainTitle, 'FontWeight', 'bold', 'FontSize', 14);
else
    % For older MATLAB versions, use a workaround
    % Add an invisible subplot to hold the main title
    subplot(2,1,1);
    title('');
    subplot(2,1,2);
    title('');
    annotation('textbox', [0.5, 0.95, 0, 0], ...
               'String', mainTitle, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'FontWeight', 'bold', ...
               'FontSize', 14);
end
end

end
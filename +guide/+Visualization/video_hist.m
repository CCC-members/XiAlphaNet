clc;
clear all;

% Paths
DataPath = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';


% Subfolders within the main folder
subFolders = {'Control'};

% Load parameters
load('Data\Model_Parameters\parameters.mat');
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Na = parameters.Dimensions.Na;
Nw = parameters.Dimensions.Nw;

%load('Data\RinvT\RinvT_G.mat', 'RinvT');

All_Data = {};  % Initialize a cell array to store data

% Loop through the subfolders and process the data
for k = 1:length(subFolders)
   % folderPath_D = fullfile(DataPath, subFolders{k});
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
   % folderPath_Age = fullfile(dataAge, subFolders{k});

    % Get a list of all .mat files in the subfolder
   % matFiles_D = dir(fullfile(folderPath_D, '*.mat'));
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
   % matFiles_Age = dir(fullfile(folderPath_Age, '*.mat'));
   
    for j = 1:length(matFiles_X)
        %filePath_D = fullfile(folderPath_D, matFiles_D(j).name);
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        %filePath_Age = fullfile(folderPath_Age, matFiles_Age(j).name);
        
        % Load the .mat files
        %data_D = load(filePath_D);
        data_X = load(filePath_X);
        %age = load(filePath_Age);
        
        % Extract and store the data
        All_Data{1,j} = data_X.x(1:end-1);
        All_Data{2,j} = data_X.x(end);
    end
end
% Initialize age partitions with interpolation
age_partitions = linspace(0, 100, 100);  % Adjust the age range as needed
interp_factor = 5;  % Interpolation factor to slow down and smooth the video
interpolated_ages = linspace(min(age_partitions), max(age_partitions), length(age_partitions) * interp_factor);

% Prepare video writer
videoFileName = 'density_video_smooth_with_full_histogram_slow.mp4';
video = VideoWriter(videoFileName, 'MPEG-4');
open(video);

% Create a single figure for plotting
figureHandle = figure;

% Percentile limits initialization
p_low = 1;  % Lower percentile for range limit
p_high = 99;  % Upper percentile for range limit
log_a_all_data = [];
log_e_all_data = [];

% First pass to collect log_a_all and log_e_all data for tighter range calculation
for p = 1:length(age_partitions)-1
    a_all = [];
    e_all = [];
    
    % Loop through all data
    for j = 1:length(All_Data(1,:))
        if All_Data{2,j} >= age_partitions(p) && All_Data{2,j} < age_partitions(p+1)
            % Extract e, a, and s2 values
            [e, a, s2] = x2v(All_Data{1,j});
            a_all = [a_all; a(:,4)];  % Store the 4th column of 'a' (Alpha amplitude)
            e_all = [e_all; a(:,1)];  % Store the 1st column of 'a' (Xi amplitude)
        end
    end
    
    if isempty(a_all) || isempty(e_all)
        continue;
    end

    % Filter and log-transform valid data points
    valid_indices = a_all > 7;
    log_a_all = log(a_all(valid_indices));
    log_e_all = log(e_all(valid_indices));

    if isempty(log_a_all) || isempty(log_e_all)
        continue;
    end

    % Collect data for percentile calculation
    log_a_all_data = [log_a_all_data; log_a_all];
    log_e_all_data = [log_e_all_data; log_e_all];
end

% Calculate tighter range using percentiles
lower_log_a = prctile(log_a_all_data, p_low);
upper_log_a = prctile(log_a_all_data, p_high);
lower_log_e = prctile(log_e_all_data, p_low);
upper_log_e = prctile(log_e_all_data, p_high);

% Initialize variables for tracking the initial point
initial_point_x = [];
initial_point_y = [];

% Smoothing window size for moving average
smooth_window_size = 5;  % Adjust as needed

% Main loop for generating the video with interpolated ages and fixed axis limits
for interp_age = interpolated_ages
    a_all = [];
    e_all = [];
    
    % Loop through all data
    for j = 1:length(All_Data(1,:))
        if All_Data{2,j} >= interp_age && All_Data{2,j} < interp_age + 1  % Interpolate between age partitions
            % Extract e, a, and s2 values
            [e, a, s2] = x2v(All_Data{1,j});
            a_all = [a_all; a(:,4)];  % Store the 4th column of 'a' (Alpha amplitude)
            e_all = [e_all; a(:,1)];  % Store the 1st column of 'a' (Xi amplitude)
        end
    end
    
    if isempty(a_all) || isempty(e_all)
        continue;
    end

    % Filter and log-transform valid data points
    valid_indices = a_all > 7;
    log_a_all = log(a_all(valid_indices));
    log_e_all = log(e_all(valid_indices));

    if isempty(log_a_all) || isempty(log_e_all)
        continue;
    end

    % Smooth the data using a moving average to reduce fluctuations
    if length(log_a_all) > smooth_window_size
        log_a_all = movmean(log_a_all, smooth_window_size);
        log_e_all = movmean(log_e_all, smooth_window_size);
    end

    % Create histogram for log_a_all and log_e_all
    clf(figureHandle);  % Clear the figure before each plot
    
    % Create 2D histogram with normalized count density and full coverage
    hist_obj = histogram2(log_a_all, log_e_all, 50, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'Normalization', 'countdensity');
    
    % Use the percentile-based axis limits for a tighter range
    xlim([lower_log_a, upper_log_a]);
    ylim([lower_log_e, upper_log_e]);
    
    % Set the tick marks and labels for the original scale
    original_a_ticks = exp(linspace(lower_log_a, upper_log_a, 10));
    original_e_ticks = exp(linspace(lower_log_e, upper_log_e, 10));

    % Set the tick marks in log scale
    set(gca, 'XTick', log(original_a_ticks));
    set(gca, 'YTick', log(original_e_ticks));

    % Set the tick labels in the original scale
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.2f', x), original_a_ticks, 'UniformOutput', false));
    set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.2f', x), original_e_ticks, 'UniformOutput', false));

    % Find the bin with the highest count (most probable point)
    [~, max_bin_idx] = max(hist_obj.Values(:));  % Index of the max bin
    [bin_x_idx, bin_y_idx] = ind2sub(size(hist_obj.Values), max_bin_idx);  % Convert linear index to bin indices

    % Get the center of the most probable bin
    most_probable_x = mean([hist_obj.XBinEdges(bin_x_idx), hist_obj.XBinEdges(bin_x_idx + 1)]);
    most_probable_y = mean([hist_obj.YBinEdges(bin_y_idx), hist_obj.YBinEdges(bin_y_idx + 1)]);
    
    % Save the first point if not already set
    if isempty(initial_point_x) || isempty(initial_point_y)
        initial_point_x = most_probable_x;
        initial_point_y = most_probable_y;
    end
    
    % Plot a red marker at the most probable point
    hold on;
    plot(most_probable_x, most_probable_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Plot the initial point with a different color (blue)
    plot(initial_point_x, initial_point_y, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    hold off;

    % Add labels and title with the interpolated age value
    xlabel('Alpha Peak Freq');
    ylabel('Alpha Amplitude');
    title(sprintf('Joint Histogram for Age %.2f', interp_age));  % Show only the interpolated age
    
    % Capture the frame for the video
    frame = getframe(figureHandle);
    
    % Write the frame multiple times to slow down the video
    for repeat_frame = 1:10  % Repeat each frame 10 times to slow down the video
        writeVideo(video, frame);
    end
end

% Close the video writer
close(video);

disp('Video generation complete.');



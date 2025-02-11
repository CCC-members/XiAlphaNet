function plot_peak_freq_vs_age_vs_region(region)

% This extracts the log spectrum and performs regression on peak alpha frequency
vertices = [];
if isempty(region) 
    % Occipital Lobe Regions
    region = 'L_V1_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V2_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V3A_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V3B_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V3CD_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V3_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V4_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    region = 'L_V7_ROI L';
    vertices = [vertices, getRegionVertices(region)];
    
%     % Add cuneus regions for the occipital lobe
%     region = 'L_V6_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_V6A_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_V8_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     
%     % Add parietal lobe regions
%     region = 'L_7m_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_7PC_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_7AL_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_7PL_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_7Pm_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_7Am_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_PGi_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_PGp_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'L_PGs_ROI L';
%     vertices = [vertices, getRegionVertices(region)];
    
%     % Repeat for the right hemisphere
%     region = 'R_V1_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V2_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V3A_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V3B_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V3CD_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V3_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V4_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V7_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     
%     % Add cuneus regions for the occipital lobe (right hemisphere)
%     region = 'R_V6_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V6A_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_V8_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     
%     % Add parietal lobe regions (right hemisphere)
%     region = 'R_7m_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_7PC_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_7AL_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_7PL_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_7Pm_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_7Am_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_PGi_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_PGp_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
%     region = 'R_PGs_ROI R';
%     vertices = [vertices, getRegionVertices(region)];
end

% Define paths
%DataPath = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';
%dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control'};

% Load necessary parameters
load('Data\Model_Parameters\parameters.mat');
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Na = parameters.Dimensions.Na;
Nw = parameters.Dimensions.Nw;

%load('Data\RinvT\RinvT_G.mat', 'RinvT');

% Initialize cell array
All_Data = {};

% Process data
for k = 1:length(subFolders)
    % Full path to the subfolder
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));

    % Loop through each .mat file in the subfolder
    index=1;
    for j =randperm(length(matFiles_X),1000)
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        data = load(filePath_X);

        [e,a,s2] = x2v(data.x(1:end-1));
        J = a(vertices,4);

        % Exclude zero values and NaNs
        J = J(J > 0 & ~isnan(J));

        % Store in All_Data
        All_Data{1,index} = J; % Peak alpha frequencies
        All_Data{2,index} = data.x(end); % Age of the subject
        index = index + 1;
    end
end

% Perform a robust measure of central tendency
% For example, use the median or a trimmed mean to avoid the influence of outliers
average_peak_alpha = cellfun(@(x) median(x), All_Data(1,:));
ages = cell2mat(All_Data(2,:));

% Remove NaN values from ages and the corresponding entries in average_peak_alpha
valid_indices = ~isnan(ages);
ages_clean = ages(valid_indices);
average_peak_alpha_clean = average_peak_alpha(valid_indices);

% Perform Gaussian Process Regression (GPR)
gprMdl = fitrgp(ages_clean', average_peak_alpha_clean', 'KernelFunction', 'squaredexponential', 'Standardize', true);

% Predict the values over a grid of age values for smoother plotting
age_range = linspace(min(ages_clean), max(ages_clean), 100)';
[ypred, ysd] = predict(gprMdl, age_range);

% Plot the original full range results
figure;
scatter(ages_clean, average_peak_alpha_clean, 'b', 'filled');
hold on;
plot(age_range, ypred, 'r-', 'LineWidth', 2);
plot(age_range, ypred + 2*ysd, 'r--', 'LineWidth', 1); % 95% confidence interval
plot(age_range, ypred - 2*ysd, 'r--', 'LineWidth', 1); % 95% confidence interval
xlabel('Age');
ylabel('Median Peak Alpha Frequency');
title('Gaussian Process Regression of Median Peak Alpha Frequency against Age');
legend({'Data', 'GPR Mean Prediction', '95% Confidence Interval'}, 'Location', 'best');
grid on;

% Explicitly set the axis limits if necessary
xlim([min(ages_clean)-5, max(ages_clean)+5]);
ylim([min(average_peak_alpha_clean)-0.5, max(average_peak_alpha_clean)+0.5]);

hold off;


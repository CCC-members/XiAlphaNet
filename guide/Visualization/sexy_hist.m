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

% Convert each data.x to [e, a, s2] and store a(:,1) and a(:,4)
a_all = [];
e_all = [];
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    a_all = [a_all; a(:,4)];  % Store the first column of a in a_all
    e_all = [e_all; a(:,1)];  % Store the fourth column of a in e_all
end

% Filter data for a_all > 7
log_a_all = log(a_all(a_all > 7));
log_e_all = log(e_all(a_all > 7));

% Create a grid for the density estimation
[xGrid, yGrid] = meshgrid(linspace(min(log_a_all), max(log_a_all), 100), ...
                          linspace(min(log_e_all), max(log_e_all), 100));

% Estimate the density
density = ksdensity([log_a_all, log_e_all], [xGrid(:), yGrid(:)]);

% Reshape the density to match the grid
density = reshape(density, size(xGrid));

% Normalize the density to create a probability distribution
density = density / sum(density(:));

% Plot the density as a smooth map
contourf(xGrid, yGrid, density, 'LineColor', 'none');
colorbar;

% Set the tick marks and labels for the original scale
original_a_ticks = exp(linspace(min(log_a_all), max(log_a_all), 10));
original_e_ticks = exp(linspace(min(log_e_all), max(log_e_all), 10));

% Set the tick marks in log scale
set(gca, 'XTick', log(original_a_ticks));
set(gca, 'YTick', log(original_e_ticks));

% Set the tick labels in the original scale
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.2f', x), original_a_ticks, 'UniformOutput', false));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.2f', x), original_e_ticks, 'UniformOutput', false));

xlabel('Alpha Amplitude');
ylabel('Xi Amplitud');
title('Joint PDF');

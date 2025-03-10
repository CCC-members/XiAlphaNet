modelParametersPath = 'Data\Model_Parameters';
dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control'};

% Load parameters
load('Data\Model_Parameters\parameters.mat');
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Na = parameters.Dimensions.Na;
Nw = parameters.Dimensions.Nw;

% Load sex vector
load('Data\Age\sex.mat');  % Assuming the sex data is stored in 'sex.mat' file
% The vector 'sex' should be binary (1 for male, 0 for female)

% Initialize data storage
allData_e = [];
allData_a = [];
allData_apf = [];
sexLabels = [];

% Loop over groups (but now we only care about collecting data and sex labels)
for k = 1:length(subFolders)
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));

    % Loop through each .mat file in the subfolder
    for j = 1:length(matFiles_X)
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        
        % Load the .mat file
        data = load(filePath_X);
        
        % Assuming 'x' is a vector or matrix in the loaded file
        if isfield(data, 'x')
            [e, a, s2] = x2v(data.x(:));  % Extract features
            allData_e = [allData_e, e(:,1)];  % Concatenate 'e' data
            allData_a = [allData_a, a(:,1)];  % Concatenate 'a' data
            allData_apf = [allData_apf, a(:,4)];  % Concatenate 'apf' data
            sexLabels = [sexLabels; sex(j)]; % Sex data corresponding to the file
        end
    end
end

% Run t-SNE for each dataset and visualize by sex
datasets = {allData_e, allData_a, allData_apf};
datasetNames = {'Energy', 'Activation', 'Peak Frequency'};
colors = ['r', 'b'];  % Colors for females (0) and males (1)
markers = ['o', '+'];  % Different markers for females and males

% Create a figure with subplots
figure;
for i = 1:length(datasets)
    % Perform t-SNE
    mappedX = tsne(datasets{i}', 'NumDimensions', 2);  % Transpose data to have variables in rows
    
    % Subplot for each dataset
    subplot(1, length(datasets), i);
    gscatter(mappedX(:,1), mappedX(:,2), sexLabels, colors, markers);
    title(['t-SNE of ', datasetNames{i}, ' by Sex']);
    legend('Female', 'Male');
end

% Enhance the plot
sgtitle('t-SNE Visualization by Sex');

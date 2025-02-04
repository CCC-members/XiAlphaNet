% Setup
clear all 
clc

%% Define the path to your main folder containing the subfolders
mainFolder = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Age';

% Subfolders within the main folder
subFolders = {'Control', 'Pathological'};

% Load previously saved parameters if necessary
load('Data\Model_Parameters\parameters.mat');

%% Ensure model parameters subfolders exist and create if they do not
for s = 1:length(subFolders)
    if ~exist(fullfile(modelParametersPath, subFolders{s}), 'dir')
        mkdir(fullfile(modelParametersPath, subFolders{s}));
    end
end
sex = [];
%% Loop through each subfolder
for k = 1:length(subFolders)
    % Full path to the subfolder

    folderPath = fullfile(mainFolder, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles = dir(fullfile(folderPath, '*.mat'));
   
    %% Loop through each .mat file in the subfolder
    for j = 1:length(matFiles)
        % Full path to the .mat file

        filePath = fullfile(folderPath, matFiles(j).name);
        % Load the .mat file
        data_struct = load(filePath);
        
        %% Extract and store the required data and parameters
        %parameters.Data_Cross = data_struct.data_struct.CrossM(:,:,1:49); % Assuming CrossM is the correct field
        if data_struct.data_struct.sex == 'M'
            index_sex= 1;
        else 
            index_sex=0;
        end

        sex = [sex,index_sex];
      
%         %% 
%         a = data_struct.data_struct.age;
%         %% Save the computed x to the corresponding group folder in Model_Parameters
%         saveFilePath = fullfile(modelParametersPath, subFolders{k}, sprintf('age_%d.mat', j));
%         save(saveFilePath, 'a');
    end
end

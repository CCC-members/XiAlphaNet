% Setup
clear all 
clc

%% Define the path to your main folder containing the subfolders
mainFolder = 'Data/Scalp_Density_Matrix';
modelParametersPath = 'Data/Model_Parameters';

% Subfolders within the main folder
subFolders = {'Control'}; 

% Load previously saved parameters
load('Data/Model_Parameters/parameters.mat');

%% Ensure model parameters subfolders exist and create if they do not
for s = 1:length(subFolders)
    if ~exist(fullfile(modelParametersPath, subFolders{s}), 'dir')
        mkdir(fullfile(modelParametersPath, subFolders{s}));
    end
end
lambda_opt_list = [];
%% Loop through each subfolder
for k = 1:length(subFolders)
    % Full path to the subfolder

    folderPath = fullfile(mainFolder, subFolders{k});
   % folderPath_x = fullfile(modelParametersPath, subFolders{k});

    % Get a list of all .mat files in the subfolder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    %matFiles_x = dir(fullfile(folderPath_x, '*.mat'));
    
    %% Loop through each .mat file in the subfolder
    perml= randperm(length(matFiles),length(matFiles));
    for j = perml;%length(matFiles)
        % Full path to the .mat file

        filePath = fullfile(folderPath, matFiles(j).name);
        %filePath_x = fullfile(folderPath_x, matFiles_x(j).name);

        % Load the .mat file
        Nw = 47;
        data_struct = load(filePath);
        freq = data_struct.data_struct.freqrange(1:Nw);

        % Cross
        Cross = data_struct.data_struct.CrossM(:,:,1:Nw);
        Cross = aveReference(Cross);
       
        %% Extract and store the required data and parameters
        parameters1  = parameters;
        parameters1.Data.Cross = Cross; 
        parameters1.Data.freq = freq;
        parameters1.Age = data_struct.data_struct.age;
        if ischar(data_struct.data_struct.age) || isstring(data_struct.data_struct.age)
            % Convert the string to a number
            age = str2double(data_struct.data_struct.age);
        else
            % If it's already a number, just use it directly
            age = data_struct.data_struct.age;
        end
        
        %% 
        
        tic
        [d] = Xi_ALphaNET_delays(parameters1);
        toc
        d.Age=age;
        %% Save the computed x to the corresponding group folder in Model_Parameters
        saveFilePath = fullfile(modelParametersPath,'Delays&Conn', sprintf('d_opt_%d.mat', j));
        save(saveFilePath, 'd');
    end
end




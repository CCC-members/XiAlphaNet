% Setup
clear all 
clc

%% Define the path to your main folder containing the subfolders
mainFolder = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';

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

%% Loop through each subfolder
for k = 1:length(subFolders)
    % Full path to the subfolder

    folderPath = fullfile(mainFolder, subFolders{k});
   % folderPath_x = fullfile(modelParametersPath, subFolders{k});

    % Get a list of all .mat files in the subfolder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    %matFiles_x = dir(fullfile(folderPath_x, '*.mat'));
    
    %% Loop through each .mat file in the subfolder
    for j = 1:length(matFiles)
        % Full path to the .mat file

        filePath = fullfile(folderPath, matFiles(j).name);
        %filePath_x = fullfile(folderPath_x, matFiles_x(j).name);

        % Load the .mat file
        data_struct = load(filePath);
        
        %% Extract and store the required data and parameters
        parameters.Data.Cross = data_struct.data_struct.CrossM(:,:,1:49); % Assuming CrossM is the correct field
        parameters.Age = data_struct.data_struct.age;
        
        %% 
        % 56022
        %k_min = findMinimumK(parameters,11,20); % Estimate the number of frequencies that approximate the f and dF with a relative error less than 5%
        parameters.Stochastic.stoch = 1;
        parameters.Stochastic.Nsfreq = 1;
        parameters.Stochastic.Niter = 1;
        parameters = sample_frequencies(parameters);
        parameters.Threshold = activation_threshold(parameters);
        parameters.Lipschitz = 400; % estimateLipschitzConstant(parameters, 1, 60);
        lambda_space = lambda_regspace(parameters);
        disp('-->> Initializing Stochastic FISTA global optimizer...')
        [lambda_opt] = bayesianOptSearch(lambda_space, parameters);
        parameters.Stochastic.stoch = 0;  
        [x_opt, hist] = stoch_fista_global(lambda_opt, parameters);
        x = x_opt.Solution;
        
        %% Save the computed x to the corresponding group folder in Model_Parameters
        saveFilePath = fullfile(modelParametersPath, subFolders{k}, sprintf('x_opt_%d.mat', j));
        save(saveFilePath, 'x');
    end
end



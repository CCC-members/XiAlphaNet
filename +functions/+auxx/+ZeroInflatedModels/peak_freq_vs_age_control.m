% This extracts the log spectrum and performs regression on peak alpha frequency
clc;
clear all;

% Define paths
DataPath = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';
dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control','Pathological'};

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
    folderPath_D = fullfile(DataPath, subFolders{k});
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    folderPath_Age = fullfile(dataAge, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles_D = dir(fullfile(folderPath_D, '*.mat'));
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
    matFiles_Age = dir(fullfile(folderPath_Age, '*.mat'));
   
    % Loop through each .mat file in the subfolder
    for j = 1:length(matFiles_X)
        %%
        filePath_D = fullfile(folderPath_D, matFiles_D(j).name);
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        filePath_Age = fullfile(folderPath_Age, matFiles_Age(j).name);
        
        % Load data
        load(filePath_D);
        c = data_struct.CrossM;
        
        data = load(filePath_X);
        age = load(filePath_Age);
        
        [e,a,s2] = x2v(data.x);
        J = a(:,4);
        
        % Select relevant voxels (example indices given, adjust as needed)
        r1  = 1;
        r2 =  3900;
        delta = 1;
        J1 = J(r1:(r1+delta));
        J2 = J(r2:(r2+delta));
        J = [J1; J2];
        
        % Store in All_Data
        All_Data{1,j} = J; % Peak alpha frequencies
        All_Data{2,j} = age.a; % Age of the subject
    end
end

% Perform the averaging of peak alpha frequencies across voxels for each subject
average_peak_alpha = cellfun(@mean, All_Data(1,:));
ages = cell2mat(All_Data(2,:));

% Perform Gaussian Process regression
gprMdl = fitrgp(ages', average_peak_alpha', 'KernelFunction', 'squaredexponential');

% Predict the values over a grid of age values for smoother plotting
age_range = linspace(min(ages), max(ages), 100)';
[ypred, ysd] = predict(gprMdl, age_range);
age_range=age_range;
% Plot the results
figure;
scatter(ages, average_peak_alpha, 'b', 'filled');
hold on;
plot(age_range, ypred, 'r-', 'LineWidth', 2);
plot(age_range, ypred + 2*ysd, 'r--', 'LineWidth', 1); % 95% confidence interval
plot(age_range, ypred - 2*ysd, 'r--', 'LineWidth', 1); % 95% confidence interval
xlabel('Age');
ylabel('Average Peak Alpha Frequency');
title('Gaussian Process Regression of Average Peak Alpha Frequency against Age');
legend({'Data', 'GP Mean Prediction', '95% Confidence Interval'}, 'Location', 'best');
grid on;
hold off;

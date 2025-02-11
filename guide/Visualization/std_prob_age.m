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
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
   
    for j = 1:length(matFiles_X)
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        data_X = load(filePath_X);
        
        % Extract and store the data
        All_Data{1,j} = data_X.x(1:end-1);
        All_Data{2,j} = data_X.x(end);
    end
end
% Initialize storage for Peak Alpha Frequency (a(:,4)), Amplitude of the Alpha (a(:,1)), and Amplitude of Xi
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);         % Store Amplitude of the Xi
end

% Calculate sparse index (optional, can be used for any of the measures)
sparse_index_PAF = sum(PAF_all(:)==0)/length(PAF_all(:))*100;
sparse_index_AlphaAmp = sum(AlphaAmp_all(:)==0)/length(AlphaAmp_all(:))*100;
sparse_index_XiAmp = sum(XiAmp_all(:)==0)/length(XiAmp_all(:))*100;

% Define threshold for each measure
threshold_PAF = 8;%prctile(PAF_all(:), 80);
threshold_AlphaAmp =0.8;%prctile(AlphaAmp_all(:), 80);
threshold_XiAmp =0.8;%prctile(XiAmp_all(:), 80);

% Threshold and binary mask for each measure
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4) > threshold_PAF;         % Apply threshold for PAF
    AlphaAmp_all(:,j) = a(:,1) > threshold_AlphaAmp; % Apply threshold for Alpha Amplitude
    XiAmp_all(:,j) = e(:,1) > threshold_XiAmp;     % Apply threshold for Xi Amplitude
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages

% Define age intervals
age_intervals = linspace(0, 100, 6);  % 5 intervals, 6 edges

% Initialize storage for the averaged values
PAF_avg_intervals = zeros(size(PAF_all, 1), length(age_intervals)-1);
AlphaAmp_avg_intervals = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
XiAmp_avg_intervals = zeros(size(XiAmp_all, 1), length(age_intervals)-1);


% Standard deviation values for each age interval
for i = 1:length(age_intervals)-1
    age_mask = (ages >= age_intervals(i)) & (ages < age_intervals(i+1));
    
    % For PAF standard deviation
    PAF_group = PAF_all(:, age_mask);
    PAF_std_intervals(:, i) = std(PAF_group, 0, 2);  % Compute std across subjects in the group
    
    % For Alpha Amplitude standard deviation
    AlphaAmp_group = AlphaAmp_all(:, age_mask);
    AlphaAmp_std_intervals(:, i) = std(AlphaAmp_group, 0, 2);  % Compute std across subjects in the group
    
    % For Xi Amplitude standard deviation
    XiAmp_group = XiAmp_all(:, age_mask);
    XiAmp_std_intervals(:, i) = std(XiAmp_group, 0, 2);  % Compute std across subjects in the group
end


% Plotting for Peak Alpha Frequency (PAF) in the first row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*XiAmp_std_intervals(:,i) ;
    J = J_age_interval;
   % subplot(3, 5, i);  % 3 rows, 5 columns, first row for PAF
    %esi_plot(gca, J_age_interval, [0, max(PAF_avg_intervals(:))]);
    %colormap("hot");
    esi_plot_single
    title(sprintf('Alpha Peak Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end


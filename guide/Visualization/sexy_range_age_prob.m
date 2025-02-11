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
x_avg=zeros(7*Nr+1,1899);
% Loop through the subfolders and process the data
for k = 1:length(subFolders)
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
   
    for j = 1:length(matFiles_X)
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        data_X = load(filePath_X);
        
        % Extract and store the data
        %x_avg(:,j) =  data_X.x.Solution;
        All_Data{1,j} = data_X.x.Solution;
        All_Data{2,j} = data_X.x.Age;
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

% Average values for each age interval
for i = 1:length(age_intervals)-1
    age_mask = (ages >= age_intervals(i)) & (ages < age_intervals(i+1));
    
    % For PAF
    f_va_PAF = sum(PAF_all(:, age_mask), 2);
    nf_va_PAF = sum(age_mask);  % Number of subjects in this age group
    PAF_avg_intervals(:, i) = f_va_PAF / nf_va_PAF;
    
    % For Alpha Amplitude
    f_va_AlphaAmp = sum(AlphaAmp_all(:, age_mask), 2);
    nf_va_AlphaAmp = sum(age_mask);  % Number of subjects in this age group
    AlphaAmp_avg_intervals(:, i) = f_va_AlphaAmp / nf_va_AlphaAmp;
    
    % For Xi Amplitude
    f_va_XiAmp = sum(XiAmp_all(:, age_mask), 2);
    nf_va_XiAmp = sum(age_mask);  % Number of subjects in this age group
    XiAmp_avg_intervals(:, i) = f_va_XiAmp / nf_va_XiAmp;
end

% Create a single large figure for all subplots
%figure('Position', [100, 100, 1500, 900]);  % Adjust figure size for better visibility

% Plotting for Peak Alpha Frequency (PAF) in the first row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*PAF_avg_intervals(:,i) ;
    J = J_age_interval;
   % subplot(3, 5, i);  % 3 rows, 5 columns, first row for PAF
    %esi_plot(gca, J_age_interval, [0, max(PAF_avg_intervals(:))]);
    %colormap("hot");
    esi_plot_single
    title(sprintf('Alpha Peak Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Plotting for Alpha Amplitude in the second row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*AlphaAmp_avg_intervals(:,i) ;
    
    %subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
    %esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
    %colormap("hot");
    J = J_age_interval;
    esi_plot_single
    title(sprintf('Alpha Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Plotting for Xi Amplitude in the third row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*XiAmp_avg_intervals(:,i) ;
    
    %subplot(3, 5, i + 10);  % 3 rows, 5 columns, third row for Xi Amplitude
    %esi_plot(gca, J_age_interval, [0, max(XiAmp_avg_intervals(:))]);
    %colormap("hot");
    J = J_age_interval;
    esi_plot_single
    title(sprintf('Xi Amplitud vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Adjust layout and appearance
sgtitle('Lifespan Mapping of Probability Distribution');  % Overall title

% Create a single large figure for all subplots
figure('Position', [100, 100, 1500, 900]);  % Adjust figure size for better visibility

% % Plotting for Peak Alpha Frequency (PAF) in the first row
% for i = 1:length(age_intervals)-1
%     J_age_interval = 1.*(PAF_avg_intervals(:,i) > 0.90);
%     J = J_age_interval;
%    % subplot(3, 5, i);  % 3 rows, 5 columns, first row for PAF
%    % esi_plot(gca, J_age_interval, [0, max(PAF_avg_intervals(:))]);
%    esi_plot_single
%     %colormap("hot");
%     title(sprintf('Alpha Peak Freq vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
% end

% Plotting for Alpha Amplitude in the second row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*AlphaAmp_avg_intervals(:,i) > 0.99;
    J = J_age_interval;
    %subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
    %esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
    esi_plot_single
    title(sprintf('Alpha Peak Freq vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
    %colormap("hot");
end

% Plotting for Xi Amplitude in the third row
for i = 1:length(age_intervals)-1
    J_age_interval = 1.*(XiAmp_avg_intervals(:,i) > 0.95);
    J = J_age_interval;
    %subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
    %esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
    esi_plot_single
    title(sprintf('Xi Amplitude vs Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
end

% Adjust layout and appearance
sgtitle('Lifespan Mapping of High-Probability Voxel Activity (p>0.5)');  % Overall title

% Create a single large figure for all subplots
figure('Position', [100, 100, 1500, 900]);  % Adjust figure size for better visibility

% % Plotting for Peak Alpha Frequency (PAF) in the first row
% for i = 1:length(age_intervals)-1
%     J_age_interval = 1.*(PAF_avg_intervals(:,i) > 0.90);
%     
%     subplot(3, 5, i);  % 3 rows, 5 columns, first row for PAF
%     esi_plot(gca, J_age_interval, [0, max(PAF_avg_intervals(:))]);
%     colormap("hot");
%     title(sprintf('Age %.2f - %.2f', age_intervals(i), age_intervals(i+1)));
% end
% 
% % Plotting for Alpha Amplitude in the second row
% for i = 1:length(age_intervals)-1
%     J_age_interval = 1.*(AlphaAmp_avg_intervals(:,i) > 0.90);
%     
%     subplot(3, 5, i + 5);  % 3 rows, 5 columns, second row for Alpha Amplitude
%     esi_plot(gca, J_age_interval, [0, max(AlphaAmp_avg_intervals(:))]);
%     colormap("hot");
% end
% 
% % Plotting for Xi Amplitude in the third row
% for i = 1:length(age_intervals)-1
%     J_age_interval = 1.*(XiAmp_avg_intervals(:,i) > 0.90);
%     
%     subplot(3, 5, i + 10);  % 3 rows, 5 columns, third row for Xi Amplitude
%     esi_plot(gca, J_age_interval, [0, max(XiAmp_avg_intervals(:))]);
%     colormap("hot");
% end

% Adjust layout and appearance
%sgtitle('Lifespan Mapping of High-Probability Voxel Activity (p>0.9)');  % Overall title

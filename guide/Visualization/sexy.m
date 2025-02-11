clc;
clear all;

% Paths
DataPath = 'Data/Scalp_Density_Matrix';
modelParametersPath = 'Data/Model_Parameters';
%dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control'};

% Load parameters
load('Data/Model_Parameters/parameters.mat');
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
Nw = parameters.Dimensions.Nw;

%load('Data/RinvT/RinvT_G.mat', 'RinvT');

All_Data = {};  % Initialize a cell array to store data

% Loop through the subfolders and process the data
for k = 1:length(subFolders)
    %folderPath_D = fullfile(DataPath, subFolders{k});
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    %folderPath_Age = fullfile(dataAge, subFolders{k});

    % Get a list of all .mat files in the subfolder
    %matFiles_D = dir(fullfile(folderPath_D, '*.mat'));
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
    %matFiles_Age = dir(fullfile(folderPath_Age, '*.mat'));
   
    for j = 1:length(matFiles_X)
        %filePath_D = fullfile(folderPath_D, matFiles_D(j).name);
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        %filePath_Age = fullfile(folderPath_Age, matFiles_Age(j).name);
        
        % Load the .mat files
        %data_D = load(filePath_D);
        data_X = load(filePath_X);
       % age = load(filePath_Age);
        
        % Extract and store the data
        All_Data{1,j} = data_X.x.Solution;
        All_Data{2,j} = data_X.x.Age;
    end
end

% Initialize storage for a(:,4) for all subjects
a_all = [];

% Convert each data.x to [e, a, s2] and store a(:,4)
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    a_all(:,j) = a(:,4);  % Store the fourth column of a in a_all
end
threshold = 8;%prctile(a_all(:),95);
%
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    a_all(:,j) = e(:,1)>threshold;  % Store the fourth column of a in a_all
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
unique_ages = unique(ages);

% Initialize storage for the averaged values
a_avg = zeros(size(a_all, 1), length(unique_ages));

% Average a(:,4) values for the same age
for i = 1:length(unique_ages)
    age_mask = (ages == unique_ages(i));
    % Sum for each voxel the number of occurrences across the subjects with the same age
    f_va = sum(a_all(:, age_mask), 2);
    % Normalize by the number of subjects in this age group
    nf_va = sum(age_mask);  % Number of subjects in this age group
    a_avg(:, i) = f_va / nf_va;
end

% Apply a moving average to smooth the transitions
windowSize = 5;  % Define the size of the moving average window (adjust as needed)
a_smooth = movmean(a_avg, windowSize, 2);

% Determine global color limits across all ages
global_min = min(a_smooth(:));
global_max = max(a_smooth(:));
colorLimits = [0, global_max];

% Initialize Video Writer
videoFilename = 'J_age_video.avi';
v = VideoWriter(videoFilename, 'Uncompressed AVI');
v.FrameRate = 10;
open(v);

% Create the video by plotting J(age) = a_smooth(:,i) for each unique age
figure;
ax = axes;  % Create an axis for plotting

for i = 1:length(unique_ages)
    J_age = a_smooth(:,i);  % Get the smoothed a(:,4) at this age
    
    % Use the esi_plot function for plotting with fixed color limits
    esi_plot(ax, J_age, colorLimits);
    colormap("hot");
    title(ax, sprintf('Age: %.2f', unique_ages(i)));
    
    % Capture the frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Finalize the Video
close(v);

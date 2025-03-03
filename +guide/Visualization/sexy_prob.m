clc;
clear all;

% Paths
DataPath = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';

% Subfolders within the main folder
subFolders = {'Control'};

% Load parameters
load(fullfile(modelParametersPath, 'parameters.mat'));
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Na = parameters.Dimensions.Na;
Nw = parameters.Dimensions.Nw;

load('Data\RinvT\RinvT_G.mat', 'RinvT');

All_Data = {};  % Initialize a cell array to store data

% Loop through the subfolders and process the data
for k = 1:length(subFolders)
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
   
    for j = 1:length(matFiles_X)
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        
        % Load the .mat files
        data_X = load(filePath_X);
        
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
    a_all(:,j) = a(:,1);  % Store the fourth column of a in a_all
    e_all(:,j) = e(:,1);  % Store the fourth column of a in a_all
    w_all(:,j) = a(:,4);  % Store the fourth column of a in a_all
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
unique_ages = unique(ages);

% Aggregate all voxel data across subjects
all_voxel_data = a_all(:);

% Compute the 95th percentile threshold for the aggregated data
threshold = prctile(all_voxel_data, 95);

% Determine active voxels for each subject based on the threshold
active_voxels = a_all > threshold;

% Initialize storage for the probability distributions across ages
prob_distributions = zeros(size(a_all, 1), length(unique_ages));

% Compute the probability distributions for each age
for i = 1:length(unique_ages)
    age_mask = (ages == unique_ages(i));  % Mask for subjects of the same age
    active_voxels_age_group = active_voxels(:, age_mask);  % Active voxels for the current age group
    
    % Compute the frequency of activation for each voxel
    activation_frequency = sum(active_voxels_age_group, 2);
    
    % Normalize the activation frequency across voxels to get the probability distribution
    prob_distributions(:,i) = activation_frequency / sum(activation_frequency);
end

% Plot the evolution of the probability distributions with age
figure;
for i = 1:length(unique_ages)
    subplot(ceil(sqrt(length(unique_ages))), ceil(sqrt(length(unique_ages))), i);
    
    num_voxels = size(prob_distributions, 1);
    sqrt_voxels = sqrt(num_voxels);
    
    if mod(sqrt_voxels, 1) == 0
        reshaped_size = [sqrt_voxels, sqrt_voxels];
        imagesc(reshape(prob_distributions(:,i), reshaped_size)); % Reshape if it's a perfect square
    else
        imagesc(prob_distributions(:,i)); % Otherwise, display as a heatmap
    end
    
    colorbar;
    title(sprintf('Age: %.2f', unique_ages(i)));
    xlabel('Voxel Index');
    ylabel('Probability');
end

% Optionally, create a video of the probability distributions over ages
videoFilename = 'Probability_Distribution_Evolution.avi';
v = VideoWriter(videoFilename, 'Uncompressed AVI');
v.FrameRate = 10;
open(v);

for i = 1:length(unique_ages)
    if mod(sqrt_voxels, 1) == 0
        imagesc(reshape(prob_distributions(:,i), reshaped_size));
    else
        imagesc(prob_distributions(:,i));
    end
    colorbar;
    title(sprintf('Age: %.2f', unique_ages(i)));
    xlabel('Voxel Index');
    ylabel('Probability');
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);

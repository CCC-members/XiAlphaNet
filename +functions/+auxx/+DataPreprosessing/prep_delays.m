%%
%Preprocessing of the Delays 
% Specify the file path

filePath = 'Data\Average_Velocity_ROISpace\MNI-HCP-MMP1.txt';

% Open the file for reading
fid = fopen(filePath, 'rt');

% Initialize an empty cell array
ids = {};

% Read each line from the file
tline = fgetl(fid);
while ischar(tline)
    ids{end+1} = tline;  % Append the line to the cell array
    tline = fgetl(fid);   % Read the next line
end

% Close the file
fclose(fid);

averageDelays_FTRACTS.parcelIDs = ids;
%%
d = real(readmatrix('Data\Average_Velocity_ROISpace\dcm_axonal_delay__median.txt', 'Delimiter', ' '));
averageDelays_FTRACTS.delays  = d;

%%
 %Assuming 'd' is your matrix with some NaN values
% Compute the median of non-NaN values in the matrix
medianValue = median(d(~isnan(d)));

% Replace all NaN values in the matrix with the computed median value
d2(isnan(d)) = medianValue;

averageDelays_FTRACTS.delaysMedian = d2;

%%
 %Assuming 'd' is your matrix with some NaN values
% Compute the mean of non-NaN values in the matrix
meanValue = mean(d(~isnan(d)));

% Replace all NaN values in the matrix with the computed mean value
d2(isnan(d)) = meanValue;

averageDelays_FTRACTS.delaysMean = d2;

% Define the directory and file name
directoryPath = 'Data\Average_Velocity_ROISpace\';
fileName = 'averageDelays_FTRACTS.mat';

% Full path construction
fullPath = fullfile(directoryPath, fileName);



% This extract the log spectrum;
clc;
clear all;
%
DataPath = 'Data\Scalp_Density_Matrix';
modelParametersPath = 'Data\Model_Parameters';
%dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control'};

% 
load('Data\Model_Parameters\parameters.mat');
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Na = parameters.Dimensions.Na;
Nw = parameters.Dimensions.Nw;
index  = 1; 
for k = 1:length(subFolders)
    % Full path to the subfolder
    folderPath_D = fullfile(DataPath, subFolders{k});
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
   
    % Get a list of all .mat files in the subfolder
    matFiles_D = dir(fullfile(folderPath_D, '*.mat'));
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
   
    %% Loop through each .mat file in the subfolder
    for j = randperm(length(matFiles_X),100)
        %%
        j
        filePath_D = fullfile(folderPath_D, matFiles_D(j).name);
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        %
        load(filePath_D);
        
        c = data_struct.CrossM;
  
        %%
        J = source_cross_dummy(c,parameters);
        log_spec = log_spectrum(J,parameters);

        %
        if ischar(data_struct.age) || isstring(data_struct.age)
            % Convert the string to a number
            age = str2double(data_struct.age);
        else
            % If it's already a number, just use it directly
            age = data_struct.age;
        end
        %

        All_Data_dummy{1,index} = log_spec;
        All_Data_dummy{2,index} =  age;
        index  =  index   + 1; 
       
    end
end
All_Data = All_Data_dummy;
% Assuming All_Data is a 2x198 cell array and parameters.Data.freq contains frequency values
ages = cell2mat(All_Data(2,:));  % Extract ages as a vector
num_channels = size(All_Data{1,1}, 1);
frequencies = parameters.Data.freq;  % Use the frequency values from parameters.Data.freq
num_frequencies = length(frequencies);

% Preallocate space for log spectrum matrix for each channel
log_spectrum_all = zeros(num_channels, num_frequencies, length(ages));

% Extract log spectra for all channels
for i = 1:length(ages)
    log_spectrum_all(:,:,i) = All_Data{1,i};  % All channels' data for a given age
end

% Now perform Linear Regression of log spectrum against age for each channel and frequency
lmModels_all = cell(num_channels, num_frequencies);  % Store the linear models for each channel and frequency
for ch = 1:num_channels
    for freq = 1:num_frequencies
        % Get the log spectrum values for this channel and frequency across all ages
        y = squeeze(log_spectrum_all(ch, freq, :));
        % Perform Linear Regression: y ~ age
        lmModels_all{ch, freq} = fitlm(ages', y);
    end
end

% Now generate the 3D plot of the regressed log spectrum for the averaged data (first plot)
% Create a meshgrid for age and frequency
age_range = linspace(min(ages), max(ages), 100);  % Generate a range of ages for plotting
[FreqGrid, AgeGrid] = meshgrid(frequencies, age_range);  % Meshgrid for frequency and age

% Compute the regressed log spectrum for the age range using the linear models for averaged data
regressed_log_spectrum_avg = zeros(numel(age_range), numel(frequencies));
for freq = 1:num_frequencies
    regressed_log_spectrum_avg(:,freq) = predict(lmModels_all{1, freq}, age_range');  % Use channel 1 for demo
end

% Plotting the first 3D plot
figure;
surf(FreqGrid, AgeGrid, regressed_log_spectrum_avg);

xlabel('Frequency (Hz)');
ylabel('Age');
zlabel('Regressed Log Spectrum');
title('3D Plot of Regressed Log Spectrum Across Age Range (Linear Regression) - Averaged Data');
shading interp;  % For smoother color transition

% Now generate the 3D plot for all channels (second plot)
figure;
hold on;
for ch = 1:num_channels
    regressed_log_spectrum_ch = zeros(numel(age_range), numel(frequencies));
    for freq = 1:num_frequencies
        regressed_log_spectrum_ch(:,freq) = predict(lmModels_all{ch, freq}, age_range');
    end
    surf(FreqGrid, AgeGrid, regressed_log_spectrum_ch);
end

xlabel('Frequency (Hz)');
ylabel('Age');
zlabel('Regressed Log Spectrum');
title('3D Plot of Regressed Log Spectrum Across Age Range (Linear Regression) - All Channels');
shading interp;  % For smoother color transition
hold off;

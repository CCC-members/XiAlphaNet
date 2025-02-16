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

load('Data\RinvT\RinvT_G.mat', 'RinvT');

for k = 1:length(subFolders)
    % Full path to the subfolder
    folderPath_D = fullfile(DataPath, subFolders{k});
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
 %   folderPath_Age = fullfile(dataAge, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles_D = dir(fullfile(folderPath_D, '*.mat'))
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
  %  matFiles_Age = dir(fullfile(folderPath_Age, '*.mat'));
   
    %% Loop through each .mat file in the subfolder
    for j = 1:length(matFiles_X)
        %%
        j
        filePath_D = fullfile(folderPath_D, matFiles_D(j).name);
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        filePath_Age = fullfile(folderPath_Age, matFiles_Age(j).name);
        % ReadCross
        
        load(filePath_D);
        c = data_struct.CrossM;
        
        % Load the .mat file
        age = load(filePath_Age);
        %

        %%
        log_spec = log_spectrum(c,parameters);

        All_Data_dummy{1,j} = log_spec;
        All_Data_dummy{2,j} =  age.a;
       
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

% Now perform Gaussian Process Regression of log spectrum against age for each channel and frequency
gprModels_all = cell(num_channels, num_frequencies);  % Store the GPR models for each channel and frequency
for ch = 1:num_channels
    for freq = 1:num_frequencies
        % Get the log spectrum values for this channel and frequency across all ages
        y = squeeze(log_spectrum_all(ch, freq, :));
        % Perform GPR: y ~ GPR(age)
        gprModels_all{ch, freq} = fitrgp(ages', y, 'BasisFunction', 'constant', 'KernelFunction', 'squaredexponential', 'FitMethod', 'exact', 'PredictMethod', 'exact');
    end
end

% Now generate the 3D plot of the regressed log spectrum for the averaged data (first plot)
% Create a meshgrid for age and frequency
age_range = linspace(min(ages), max(ages), 100);  % Generate a range of ages for plotting
[FreqGrid, AgeGrid] = meshgrid(frequencies, age_range);  % Meshgrid for frequency and age

% Compute the regressed log spectrum for the age range using the GPR models for averaged data
regressed_log_spectrum_avg = zeros(numel(age_range), numel(frequencies));
for freq = 1:num_frequencies
    regressed_log_spectrum_avg(:,freq) = predict(gprModels_all{1, freq}, age_range');  % Use channel 1 for demo
end

% Plotting the first 3D plot
figure;
surf(FreqGrid, AgeGrid, regressed_log_spectrum_avg);

xlabel('Frequency (Hz)');
ylabel('Age');
zlabel('Regressed Log Spectrum');
title('3D Plot of Regressed Log Spectrum Across Age Range (GPR) - Averaged Data');
shading interp;  % For smoother color transition

% Now generate the 3D plot for all channels (second plot)
figure;
hold on;
for ch = 1:num_channels
    regressed_log_spectrum_ch = zeros(numel(age_range), numel(frequencies));
    for freq = 1:num_frequencies
        regressed_log_spectrum_ch(:,freq) = predict(gprModels_all{ch, freq}, age_range');
    end
    surf(FreqGrid, AgeGrid, regressed_log_spectrum_ch);
end

xlabel('Frequency (Hz)');
ylabel('Age');
zlabel('Regressed Log Spectrum');
title('3D Plot of Regressed Log Spectrum Across Age Range (GPR) - All Channels');
shading interp;  % For smoother color transition
hold off;

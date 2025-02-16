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
T = parameters.Model.T;
 K=parameters.Model.K;
 R = voxel_roi_map;
 K =pinv(K);
 U_map =R*K;
 for j=1:Nw
     T_omega(:,:,j) = U_map*T(:,:,j);
 end
 parameters.Model.U= T_omega;
 index  = 1;
for k = 1:length(subFolders)
    % Full path to the subfolder
    folderPath_D = fullfile(DataPath, subFolders{k});
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
   % folderPath_Age = fullfile(dataAge, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles_D = dir(fullfile(folderPath_D, '*.mat'));
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
    %matFiles_Age = dir(fullfile(folderPath_Age, '*.mat'));
   
    %% Loop through each .mat file in the subfolder
    for j = randperm(length(matFiles_X),100)
        %%
        j
        filePath_D = fullfile(folderPath_D, matFiles_D(j).name);
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        %filePath_Age = fullfile(folderPath_Age, matFiles_Age(j).name);
        % ReadCross
        
        load(filePath_D);
        c = data_struct.CrossM;
        
        % Load the .mat file
        data = load(filePath_X);
       % age = load(filePath_Age);
        %
        xx=data.x(1:end-1);
        
        [c] = eval_source_conn(xx, parameters);
        
        log_spec = log_spectrum(c,parameters);
        All_Data{1,index}=  log_spec;
        All_Data{2,index}= data.x(end);
        index  = index + 1;
        %%
%         tic
%         
%         for i=1:Nw
%             J(:,:,i) = U_map*c(:,:,i)*U_map';
%         end
%         log_spec = log_spectrum(J,parameters);
% 
%         All_Data_dummy{1,j} = log_spec;
%         All_Data_dummy{2,j} =  age.a;
%        toc
    end
end
%All_Data = All_Data_dummy;
% % Assuming All_Data is a 2x198 cell array and parameters.Data.freq contains frequency values
% ages = cell2mat(All_Data(2,:));  % Extract ages as a vector
% num_channels = 360;%size(All_Data{1,1}, 1);
% frequencies = parameters.Data.freq;  % Use the frequency values from parameters.Data.freq
% num_frequencies = length(frequencies);
% 
% % Preallocate space for log spectrum matrix for each channel
% log_spectrum_all = zeros(num_channels, num_frequencies, length(ages));
% 
% % Extract log spectra for all channels
% for i = 1:length(ages)
%     log_spectrum_all(:,:,i) = All_Data{1,i};  % All channels' data for a given age
% end
% 
% 
% % Now perform Gaussian Process Regression of log spectrum against age for each channel and frequency
% gprModels_all = cell(num_channels, num_frequencies);  % Store the GPR models for each channel and frequency
% for ch = 1:num_channels
%     for freq = 1:num_frequencies
%         % Get the log spectrum values for this channel and frequency across all ages
%         y = squeeze(log_spectrum_all(ch, freq, :));
%         % Perform GPR: y ~ GPR(age)
%         gprModels_all{ch, freq} = fitrgp(ages', y, 'BasisFunction', 'constant', 'KernelFunction', 'squaredexponential', 'FitMethod', 'exact', 'PredictMethod', 'exact');
%     end
% end
% 
% % Now generate the 3D plot of the regressed log spectrum for the averaged data (first plot)
% % Create a meshgrid for age and frequency
% age_range = linspace(min(ages), max(ages), 100);  % Generate a range of ages for plotting
% [FreqGrid, AgeGrid] = meshgrid(frequencies, age_range);  % Meshgrid for frequency and age
% 
% % Compute the regressed log spectrum for the age range using the GPR models for averaged data
% regressed_log_spectrum_avg = zeros(numel(age_range), numel(frequencies));
% for freq = 1:num_frequencies
%     regressed_log_spectrum_avg(:,freq) = predict(gprModels_all{1, freq}, age_range');  % Use channel 1 for demo
% end
% 
% % Plotting the first 3D plot
% figure;
% surf(FreqGrid, AgeGrid, regressed_log_spectrum_avg);
% 
% xlabel('Frequency (Hz)');
% ylabel('Age');
% zlabel('Regressed Log Spectrum');
% title('3D Plot of Regressed Log Spectrum Across Age Range (GPR) - Averaged Data');
% shading interp;  % For smoother color transition


%%
% Assuming All_Data is a 2x198 cell array and parameters.Data.freq contains frequency values
ages = cell2mat(All_Data(2,:));  % Extract ages as a vector
num_channels = 360; % Set the number of channels
frequencies = parameters.Data.freq;  % Use the frequency values from parameters.Data.freq
num_frequencies = length(frequencies);

% Initialize list to keep track of valid data indices
valid_indices = [];

% Initialize an empty cell array to hold valid log spectra
valid_log_spectra = {};

% Extract log spectra for all channels
for i = 1:length(ages)
    data_i = All_Data{1,i};
    if isempty(data_i)
        fprintf('Warning: All_Data{1,%d} is empty and will be skipped.\n', i);
    else
        [rows, cols] = size(data_i);
        if rows == num_channels && cols == num_frequencies
            valid_log_spectra{end+1} = data_i;  % Store valid data
            valid_indices(end+1) = i;           % Store valid index
        else
            fprintf('Warning: All_Data{1,%d} is of size %d x %d, expected %d x %d. Skipping this entry.\n', i, rows, cols, num_channels, num_frequencies);
        end
    end
end

% Update ages to include only valid entries
valid_ages = ages(valid_indices);

% Convert valid_log_spectra to a 3D array
num_valid = length(valid_log_spectra);
log_spectrum_all = zeros(num_channels, num_frequencies, num_valid);
for i = 1:num_valid
    log_spectrum_all(:,:,i) = valid_log_spectra{i};
end

% Now perform linear regression of log spectrum against age for each channel and frequency
coeffs_all = zeros(num_channels, num_frequencies, 2);  % Store the linear model coefficients for each channel and frequency

for ch = 1:num_channels
    for freq = 1:num_frequencies
        % Get the log spectrum values for this channel and frequency across all valid ages
        y = squeeze(log_spectrum_all(ch, freq, :));
        % Perform linear regression: y ~ age
        X = [ones(length(valid_ages), 1) valid_ages'];  % Design matrix for linear regression (intercept + age)
        coeffs = X \ y;  % Solve the linear system to get the coefficients
        coeffs_all(ch, freq, :) = coeffs;  % Store the coefficients (intercept and slope)
    end
end

% Now generate the 3D plot of the regressed log spectrum for the averaged data (first plot)
% Create a meshgrid for age and frequency
age_range = linspace(min(valid_ages), max(valid_ages), 100);  % Generate a range of ages for plotting
[FreqGrid, AgeGrid] = meshgrid(frequencies, age_range);  % Meshgrid for frequency and age

% Compute the regressed log spectrum for the age range using the linear model for averaged data
regressed_log_spectrum_avg = zeros(numel(age_range), numel(frequencies));
for freq = 1:num_frequencies
    % Use the linear model for channel 1 as an example
    intercept = coeffs_all(1, freq, 1);
    slope = coeffs_all(1, freq, 2);
    regressed_log_spectrum_avg(:,freq) = intercept + slope * age_range';  % Linear model prediction
end

% Plotting the 3D plot
figure;
surf(FreqGrid, AgeGrid, regressed_log_spectrum_avg);

xlabel('Frequency (Hz)');
ylabel('Age');
zlabel('Regressed Log Spectrum');
title('3D Plot of Regressed Log Spectrum Across Age Range (Linear Fit) - Averaged Data');
shading interp;  % For smoother color transition

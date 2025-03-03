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
Nv = parameters.Dimensions.Nv;
Nw = parameters.Dimensions.Nw;

%load('Data\RinvT\RinvT_G.mat', 'RinvT');

All_Data = {};  % Initialize a cell array to store data
x_avg=zeros(7*Nv+1,1899);
% Loop through the subfolders and process the data
for k = 1:length(subFolders)
    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
   
    for j = 1:100%length(matFiles_X)
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
    PAF_all(:,j) = a(:,4).*(a(:,1)>1.8632);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1).*(a(:,1)>1.8632);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);         % Store Amplitude of the Xi
end

% % Calculate sparse index (optional, can be used for any of the measures)
% sparse_index_PAF = sum(PAF_all(:)==0)/length(PAF_all(:))*100;
% sparse_index_AlphaAmp = sum(AlphaAmp_all(:)==0)/length(AlphaAmp_all(:))*100;
% sparse_index_XiAmp = sum(XiAmp_all(:)==0)/length(XiAmp_all(:))*100;
% 
% % Define threshold for each measure
% threshold_PAF = 7;%prctile(PAF_all(:), 80);
% threshold_AlphaAmp = 1.8632;%prctile(AlphaAmp_all(:), 80);
% threshold_XiAmp =1.8632;%prctile(XiAmp_all(:), 80);
% 
% % Threshold and binary mask for each measure
% for j = 1:length(All_Data(1,:))
%     j
%     [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
%     PAF_all(:,j) = a(:,4) > threshold_PAF;         % Apply threshold for PAF
%     threshold_AlphaAmp =  1.8632;
%     AlphaAmp_all(:,j) = a(:,1) > threshold_AlphaAmp; % Apply threshold for Alpha Amplitude
%     XiAmp_all(:,j) = e(:,1) > threshold_AlphaAmp;%threshold_XiAmp;     % Apply threshold for Xi Amplitude
% end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages

% Define age intervals
age_intervals = linspace(0, 100, 6);  % 5 intervals, 6 edges


%%
% Initialize storage for conditional probabilities and expectations
PAF_cond = zeros(size(PAF_all, 1), length(age_intervals)-1);
PAF_one_minus_Pzero = zeros(size(PAF_all, 1), length(age_intervals)-1);
PAF_mean_expectation = zeros(size(PAF_all, 1), length(age_intervals)-1);

AlphaAmp_cond = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
AlphaAmp_one_minus_Pzero = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
AlphaAmp_mean_expectation = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);

XiAmp_cond = zeros(size(XiAmp_all, 1), length(age_intervals)-1);
XiAmp_one_minus_Pzero = zeros(size(XiAmp_all, 1), length(age_intervals)-1);
XiAmp_mean_expectation = zeros(size(XiAmp_all, 1), length(age_intervals)-1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zero-inflation Model Fitting Across Age Intervals %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize storage for conditional probabilities and expectations
PAF_cond              = zeros(size(PAF_all, 1), length(age_intervals)-1);
PAF_one_minus_Pzero   = zeros(size(PAF_all, 1), length(age_intervals)-1);
PAF_mean_expectation  = zeros(size(PAF_all, 1), length(age_intervals)-1);

AlphaAmp_cond         = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
AlphaAmp_one_minus_Pzero = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
AlphaAmp_mean_expectation = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);

XiAmp_cond            = zeros(size(XiAmp_all, 1), length(age_intervals)-1);
XiAmp_one_minus_Pzero = zeros(size(XiAmp_all, 1), length(age_intervals)-1);
XiAmp_mean_expectation = zeros(size(XiAmp_all, 1), length(age_intervals)-1);

for i = 1 : (length(age_intervals) - 1)
    i
    age_mask = (ages >= age_intervals(i)) & (ages < age_intervals(i+1));
    
    % Check how many subjects fall in this interval
    if sum(age_mask) == 0
        continue
    end
    
    %--------------------------------------------------------------
    % 1) Select the ages for these subjects
    %    e.g.,  if age_mask picks out M subjects,
    %    then x_in_interval is length M.
    %--------------------------------------------------------------
    x_in_interval = ages(age_mask);          % M x 1
    
    %--------------------------------------------------------------
    % 2) Select the data for these subjects:
    %    PAF_all, AlphaAmp_all, XiAmp_all are (nVoxels x nSubjects).
    %    So applying  (:, age_mask) yields
    %      (nVoxels x M)
    %--------------------------------------------------------------
    PAF_y_matrix       = PAF_all(:, age_mask);      
    AlphaAmp_y_matrix  = AlphaAmp_all(:, age_mask);
    XiAmp_y_matrix     = XiAmp_all(:, age_mask);
    
    %--------------------------------------------------------------
    % 3) Replicate the ages so that each voxel for subject j
    %    matches the correct age.  We have:
    %      nVoxels rows  x  M columns in PAF_y_matrix
    %    We want an X of that same size.
    %--------------------------------------------------------------
    nVoxels = size(PAF_y_matrix, 1);
    % Replicate each subject’s age across all nVoxels
    % x_expanded will be (nVoxels x M)
    x_expanded = repmat(x_in_interval(:).', nVoxels, 1);
    
    % Flatten these for the model:
    X_final_PAF = x_expanded(:);         % (nVoxels*M) x 1
    Y_final_PAF = PAF_y_matrix(:);       % (nVoxels*M) x 1
    
    X_final_Alpha = x_expanded(:);
    Y_final_Alpha = AlphaAmp_y_matrix(:);
    
    X_final_Xi = x_expanded(:);
    Y_final_Xi = XiAmp_y_matrix(:);
    
    %--------------------------------------------------------------
    % 4) We do NOT use random effects, so pass group = []
    %    For x_eval, we can pass the same X_final_PAF to get
    %    predictions at the same points. Alternatively, you could pass
    %    a new range for predictions if desired (e.g., 1:90).
    %--------------------------------------------------------------
    %--- Fit zero-inflated model for PAF
    model_output_PAF = fit_zero_inflation_random_effects( ...
        X_final_PAF,     ... predictor
        Y_final_PAF,     ... response
        [],              ... group = [] => no random effects
        X_final_PAF);    ... evaluate at these same points
    
    % Save results in your arrays
    %   Note: The model’s cond, zi, etc. will have length = (nVoxels*M)
    %   Usually, you might want to reshape or average across subjects...
    %   But for demonstration, we’ll do a simple reshape back to (nVoxels x M)
    condVals = model_output_PAF.predictions.cond;
    ziVals   = model_output_PAF.predictions.zi;
    
    condVals_reshaped = reshape(condVals, [nVoxels, sum(age_mask)]);
    ziVals_reshaped   = reshape(ziVals,   [nVoxels, sum(age_mask)]);
    
    % You might want to average or otherwise combine these. For example:
    PAF_cond(:, i) = mean(condVals_reshaped, 2);  % average across the M columns
    PAF_one_minus_Pzero(:, i) = mean(1 - ziVals_reshaped, 2);
    PAF_mean_expectation(:, i) = mean(condVals_reshaped .* (1 - ziVals_reshaped), 2);

    %--- Fit zero-inflated model for Alpha Amplitude
    model_output_AlphaAmp = fit_zero_inflation_random_effects( ...
        X_final_Alpha, Y_final_Alpha, [], X_final_Alpha );
    
    condValsA = model_output_AlphaAmp.predictions.cond;
    ziValsA   = model_output_AlphaAmp.predictions.zi;
    
    condValsA_reshaped = reshape(condValsA, [nVoxels, sum(age_mask)]);
    ziValsA_reshaped   = reshape(ziValsA,   [nVoxels, sum(age_mask)]);
    
    AlphaAmp_cond(:, i) = mean(condValsA_reshaped, 2);
    AlphaAmp_one_minus_Pzero(:, i) = mean(1 - ziValsA_reshaped, 2);
    AlphaAmp_mean_expectation(:, i) = mean(condValsA_reshaped .* (1 - ziValsA_reshaped), 2);

    %--- Fit zero-inflated model for Xi Amplitude
    model_output_XiAmp = fit_zero_inflation_random_effects( ...
        X_final_Xi, Y_final_Xi, [], X_final_Xi );
    
    condValsX = model_output_XiAmp.predictions.cond;
    ziValsX   = model_output_XiAmp.predictions.zi;
    
    condValsX_reshaped = reshape(condValsX, [nVoxels, sum(age_mask)]);
    ziValsX_reshaped   = reshape(ziValsX,   [nVoxels, sum(age_mask)]);
    
    XiAmp_cond(:, i) = mean(condValsX_reshaped, 2);
    XiAmp_one_minus_Pzero(:, i) = mean(1 - ziValsX_reshaped, 2);
    XiAmp_mean_expectation(:, i) = mean(condValsX_reshaped .* (1 - ziValsX_reshaped), 2);
end




%%
% Create a single figure with subplots
figure('Position', [100, 100, 1800, 1200]);  % Adjust figure size

for i = 1:length(age_intervals)-1
    % Plot for PAF
    subplot(3, length(age_intervals)-1, i);
    J = PAF_cond(:, i);
    esi_plot_single;
    title(sprintf('PAF: Cond (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
    
    subplot(3, length(age_intervals)-1, i + length(age_intervals)-1);
    J = PAF_one_minus_Pzero(:, i);
    esi_plot_single;
    title(sprintf('PAF: 1-P[Y=0] (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
    
    subplot(3, length(age_intervals)-1, i + 2*(length(age_intervals)-1));
    J = PAF_mean_expectation(:, i);
    esi_plot_single;
    title(sprintf('PAF: Mean Exp (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
end

sgtitle('Lifespan Mapping: Peak Alpha Frequency (PAF)');

% Repeat for Alpha Amplitude
figure('Position', [100, 100, 1800, 1200]);  

for i = 1:length(age_intervals)-1
    subplot(3, length(age_intervals)-1, i);
    J = AlphaAmp_cond(:, i);
    esi_plot_single;
    title(sprintf('AlphaAmp: Cond (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
    
    subplot(3, length(age_intervals)-1, i + length(age_intervals)-1);
    J = AlphaAmp_one_minus_Pzero(:, i);
    esi_plot_single;
    title(sprintf('AlphaAmp: 1-P[Y=0] (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
    
    subplot(3, length(age_intervals)-1, i + 2*(length(age_intervals)-1));
    J = AlphaAmp_mean_expectation(:, i);
    esi_plot_single;
    title(sprintf('AlphaAmp: Mean Exp (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
end

sgtitle('Lifespan Mapping: Alpha Amplitude');

% Repeat for Xi Amplitude
figure('Position', [100, 100, 1800, 1200]);  

for i = 1:length(age_intervals)-1
    subplot(3, length(age_intervals)-1, i);
    J = XiAmp_cond(:, i);
    esi_plot_single;
    title(sprintf('XiAmp: Cond (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
    
    subplot(3, length(age_intervals)-1, i + length(age_intervals)-1);
    J = XiAmp_one_minus_Pzero(:, i);
    esi_plot_single;
    title(sprintf('XiAmp: 1-P[Y=0] (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
    
    subplot(3, length(age_intervals)-1, i + 2*(length(age_intervals)-1));
    J = XiAmp_mean_expectation(:, i);
    esi_plot_single;
    title(sprintf('XiAmp: Mean Exp (%.2f-%.2f)', age_intervals(i), age_intervals(i+1)));
end

sgtitle('Lifespan Mapping: Xi Amplitude');

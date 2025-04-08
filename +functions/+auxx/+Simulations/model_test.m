%Povilas Karvelis (2024). daviolinplot - beautiful violin and raincloud plots 
% (https://github.com/povilaskarvelis/DataViz/releases/tag/v3.2.4), GitHub. Retrieved September 20, 2024.

clear all 
clc;
% Enter direction to the structural data in Results
load("/mnt/Store/Ronaldo/dev/Data/Results/structural/parameters.mat")
% Setting the parameters to preform simulation in ROI space
parameters.Model = Compact_Model;
parameters.Compact_Model = Compact_Model;
parameters.Dimensions = Dimensions;
%
clear Compact_Model
clear Model
clear Dimensions
import functions.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.CrossSpectrum.*
import functions.auxx.ModelVectorization.*
import functions.auxx.Simulations.*
import tools.*
import  functions.auxx.DataPreprosessing.*
import functions.auxx.DataPreprosessing.aveReference.*

% Setting initialiation of te simulations  
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
parameters.Dimensions.Nv = Nr;
Nv = parameters.Dimensions.Nv; % Set voxel to roi
Nw = parameters.Dimensions.Nw;
Nsim = 10;  % Number of simulated crosspectrum 
% Simulated Cross
Svv_cross=zeros(Ne,Ne,Nw,Nsim);         % Simulated observed cross-spectrum at the scalp 
Sjj_cross=zeros(Nr,Nr,Nw,Nsim);         % Simulated source cross-spectrum at ROIs
% Estimated Cross
eL_Sjj_cross= zeros(Nr,Nr,Nw,Nsim);     % Estimated source cross-spectrum using eLORETA
lcmv_Sjj_cross= zeros(Nr,Nr,Nw,Nsim);   % Estimated source cross-spectrum using LCMV
XA_Sjj_cross= zeros(Nr,Nr,Nw,Nsim);     % Estimated source cross-spectrum using Xi-AlphaNET
higgs_Sjj_cross = zeros(Nr,Nr,Nw,Nsim); % Estimated source cross-spectrum using Higgs
%
L = parameters.Compact_Model.K;
% Generate a random sample of plausible crosspectrum in the scalp 
dir_data = '/mnt/Store/Ronaldo/dev/Data/norms';
% List all the subject folders in the directory
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
selected_folders = subject_folders(randperm(length(subject_folders), Nsim));
disp("--> Simulate source and scalp cross")
N_wishart = 10;
parfor j=1:Nsim
    tic
    % Get the path for the current subject's folder and corresponding .mat file
    subject_folder = selected_folders(j).name;
    mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
   
    % Load the .mat file
    data_struct = load(mat_file_path);
    % Extract the cross-spectrum
    Svv = data_struct.data_struct.CrossM(:,:,1:Nw);
    Svv = aveReference(Svv);
    freq = data_struct.data_struct.freqrange(1:Nw);
    % Apply the mn_cross function (assuming it's a predefined function)
    Sjj = mn_cross(Svv, parameters);  % Replace parameters with actual ones
    
    % Loop to generate cross-spectra using the Wishart distribution
    for i = 1:Nw
        % Generate complex Wishart matrices
        Sjj_cross(:,:,i,j) = generate_complex_wishart(Sjj(:,:,i), N_wishart);
        
        % Apply the L matrix transformation (ensure L is defined)
        Svv_cross(:,:,i,j) = L * Sjj_cross(:,:,i,j) * L';  % Ensure L is appropriately defined
    end
    toc;
end
subject_folder = selected_folders(1).name;
 mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
 data_struct = load(mat_file_path);
 freq = data_struct.data_struct.freqrange(1:Nw);
%% 
% Init Xi-AlphaNET properties
disp("--> Estimating source cross ")

properties.model_params.nFreqs = parameters.Dimensions.Nw;
properties.model_params.BayesIter_Delay = 1;
properties.model_params.BayesIter_Reg1 = 1;
properties.model_params.BayesIter_Reg2 = 30;
properties.model_params.Nrand1 = 1;
properties.model_params.Nrand2 = 1;
properties.model_params.delay.lambda_space_cd = [[0.99,1];[0.99,1]];%[[0.4,1.6];[0.01,2]];
properties.general_params.parallel.conn_delay = 1;
properties.model_params.stoch1 = 0;
properties.model_params.stoch2 =0;
properties.model_params.tensor_field.default =0;
import app.*
import app.functions.*
import functions.*
import functions.StochasticFISTA.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.CrossSpectrum.*
import functions.auxx.ModelVectorization.*
import functions.auxx.Regularization.*
import functions.auxx.TOperator.*
import tools.*
import functions.auxx.Simulations.*
import functions.auxx.Simulations.inverse.*
import functions.auxx.Simulations.private.*

for j=1:Nsim 
    j
    tic
    %%
    data.Cross = Svv_cross(:,:,:,j);
    data.age = 25; % 
    data.freq = freq;
    [x,T] =  Xi_ALphaNET(properties,data,parameters);
    [source_act_cross] = eval_source_conn(x.Solution, data.freq,T,parameters.Model.K,parameters.Model.R,properties);
    clear x
    clear T
    XA_Sjj_cross(:,:,:,j) = source_act_cross.Cross.Full;
    %%
    parfor i=1:Nw
        i
        % Find the etimation using eLoreta and LCMV 
        source = inverse(Svv_cross(:,:,i,j),L);
        eL_Sjj_cross(:,:,i,j) = source.eloreata.Sjj;
        lcmv_Sjj_cross(:,:,i,j) = source.lcmv.Sjj;
        %higgs_Sjj_cross(:,:,i,j) = call_higgs(Svv_cross(:,:,i,j),L);
    end
    toc
end
%%



%% Read Data  and validate spectra
% Svv_cross = data_struct.Svv_cross;
% Sjj_cross = data_struct.Sjj_cross;
% eL_Sjj_cross = data_struct.eL_Sjj_cross;
% lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;
% XA_Sjj_cross = data_struct.XA_Sjj_cross;
% data_struct2.Svv_cross=Svv_cross;
% data_struct2.Sjj_cross=Sjj_cross;
% data_struct2.eL_Sjj_cross=eL_Sjj_cross;
% data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
% data_struct2.XA_Sjj_cross=XA_Sjj_cross;

% %%
% Define your parameters
Nr = size(Sjj_cross, 1); % Number of regions
Nw = Nw;%size(Sjj_cross, 3); % Number of frequencies
%Nsim = 50; % Number of simulations

% Precompute the upper triangular indices (excluding the diagonal)
idx_triu = triu(true(Nr), 1); % Only needs to be calculated once

% Preallocate memory for scores and distances based on estimated sizes
num_elements = sum(idx_triu(:)); % Number of upper triangular elements
total_size = Nsim * Nw * num_elements;

all_labels = zeros(total_size, 1);
all_eL_scores = zeros(total_size, 1);
all_lcmv_scores = zeros(total_size, 1);
all_XA_scores = zeros(total_size, 1);
%all_higgs_scores = zeros(total_size, 1);


% Preallocate distances
eL_distances_fro = zeros(Nsim * Nw, 1);
lcmv_distances_fro = zeros(Nsim * Nw, 1);
XA_distances_fro = zeros(Nsim * Nw, 1);
%higgs_distances_fro = zeros(Nsim * Nw, 1);


eL_distances_relative = zeros(Nsim * Nw, 1);
lcmv_distances_relative = zeros(Nsim * Nw, 1);
XA_distances_relative = zeros(Nsim * Nw, 1);
%higgs_distances_relative = zeros(Nsim * Nw, 1);


eL_distances_corr = zeros(Nsim * Nw, 1);
lcmv_distances_corr = zeros(Nsim * Nw, 1);
XA_distances_corr = zeros(Nsim * Nw, 1);
%higgs_distances_corr = zeros(Nsim * Nw, 1);

eL_distances_l1 = zeros(Nsim * Nw, 1); 
lcmv_distances_l1 = zeros(Nsim * Nw, 1);
XA_distances_l1 = zeros(Nsim * Nw, 1);
%higgs_distances_l1 = zeros(Nsim * Nw, 1);

% Counters for storing data efficiently
score_counter = 1;
distance_counter = 1;

for n = 1:10
    disp(n); % Display the simulation number
    for f = 1:Nw
        % Extract true and estimated covariance matrices (pre-compute abs)
        Sjj_true = abs((Sjj_cross(:, :, f, n)));
        eL_Sjj_est = abs((eL_Sjj_cross(:, :, f, n)));
        lcmv_Sjj_est = abs((lcmv_Sjj_cross(:, :, f, n)));
        XA_Sjj_est = abs((XA_Sjj_cross(:, :, f, n)));
       % higgs_Sjj_est = angle((higgs_Sjj_cross(:, :, f, n)));

        % Flatten the upper triangular parts of the matrices
        Sjj_true_vec = diag(Sjj_true);%(idx_triu));
        eL_Sjj_est_vec = diag(eL_Sjj_est);%(idx_triu));
        lcmv_Sjj_est_vec = diag(lcmv_Sjj_est);%(idx_triu));
        XA_Sjj_est_vec = diag(XA_Sjj_est);%(idx_triu));
        %higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));

        % Threshold for ground truth labels
        threshold = mean(Sjj_true_vec);
        labels = Sjj_true_vec > threshold;

        % Store labels and scores
        num_elements = numel(Sjj_true_vec);
        all_labels(score_counter:score_counter + num_elements - 1) = labels;
        all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
        all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
        all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
       % all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
        score_counter = score_counter + num_elements;

        % Compute Frobenius norm and relative error
        eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
        lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
        XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
        %higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');

        eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        %higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');

        % Compute correlation distances
        eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
        lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
        XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
        %higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
        % L1 norm
        eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
        lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
        XA_distances_l1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
        %higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
        distance_counter = distance_counter + 1;
    end
end


% Compute ROC curves using all data
[eL_X, eL_Y, ~, eL_AUC] = perfcurve(all_labels, all_eL_scores, true);
[lcmv_X, lcmv_Y, ~, lcmv_AUC] = perfcurve(all_labels, all_lcmv_scores, true);
[XA_X, XA_Y, ~, XA_AUC] = perfcurve(all_labels, all_XA_scores, true);
%[higgs_X, higgs_Y, ~, higgs_AUC] = perfcurve(all_labels, all_higgs_scores, true);




%here
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*
% Prepare data for daviolinplot of distances
methods = {'eLORETA', 'LCMV', 'XA'};%, 'higgs'};

% Ensure data are column vectors and apply the log(1 + x) transformation
fro_distances = {log(1 + eL_distances_fro(:)), log(1 + lcmv_distances_fro(:)), log(1 + XA_distances_fro(:))};%, log(1 + higgs_distances_fro(:))};
relative_distances = {log(1 + eL_distances_relative(:)), log(1 + lcmv_distances_relative(:)), log(1 + XA_distances_relative(:))};%, log(1 + higgs_distances_relative(:))};
corr_distances = {log(1+ eL_distances_corr(:)),log(1+lcmv_distances_corr(:)), log(1+XA_distances_corr(:))};%, log(1+higgs_distances_corr(:))};
l1_distances = {log(1+eL_distances_l1(:)), log(1+lcmv_distances_l1(:)),log(1+ XA_distances_l1(:))};%, log(1+ higgs_distances_l1(:))};

% Define colors for the methods (adjust these RGB values as desired)
colors = [0.2 0.6 0.8;   % Color for eLORETA
          0.8 0.4 0.2;   % Color for LCMV
          0.6 0.8 0.2;    % Color for XA
          0.5 0.5 0.5];  %Color higgs

% Plot daviolinplot for Frobenius norm distances without scatter data
% Prepare a new figure for all subplots in a 1x4 layout (excluding Correlation Distance)
figure;

% Define font sizes
titleFontSize = 18;
labelFontSize = 14;
legendFontSize = 8;

% Subplot 1: ROC Curve
subplot(1, 5, 1);
plot(eL_X, eL_Y, 'Color', colors(1,:), 'LineWidth', 2); hold on; % eLORETA
plot(lcmv_X, lcmv_Y, 'Color', colors(2,:), 'LineWidth', 2);      % LCMV
plot(XA_X, XA_Y, 'Color', colors(3,:), 'LineWidth', 2);          % XA
%plot(higgs_X, higgs_Y, 'Color', colors(4,:), 'LineWidth', 2);    % Higgs
legend(['eLORETA (AUC = ' num2str(eL_AUC, '%.2f') ')'], ...
      ['LCMV (AUC = ' num2str(lcmv_AUC, '%.2f') ')'], ...
      ['XA (AUC = ' num2str(XA_AUC, '%.2f') ')'], 'Location', 'Best', ...
      'FontSize', 10, 'FontWeight', 'bold');  % Set legend to bold and size 10
xlabel('False Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');  % Set axes labels to bold and size 14
ylabel('True Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');
title('Combined ROC', 'FontSize', 18, 'FontWeight', 'bold');  % Title to bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks and y-ticks to bold size 14
hold off;

% Subplot 2: Frobenius Norm Distance
subplot(1, 5, 2);
daviolinplot(fro_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('Frobenius Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('Frobenius Norm', 'FontSize', 18, 'FontWeight', 'bold');       % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 3: Relative Frobenius Error
subplot(1, 5, 3);
daviolinplot(relative_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('Frob Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('Frob Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 4: L1 Distance
subplot(1, 5, 4);
daviolinplot(l1_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('L1 Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('L1 Norm', 'FontSize', 18, 'FontWeight', 'bold');        % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 5: Spectral Norm Distance (L2 Distance)
subplot(1, 5, 5);
daviolinplot(corr_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('L1 Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('L1 Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Adjust the layout for better visualization
set(gcf, 'Position', [100, 100, 1800, 500]); % Adjust figure size for the 1x5 wide layout
set(gcf,'Color','w');

%% Read Data 
% Svv_cross = data_struct.Svv_cross;
% Sjj_cross = data_struct.Sjj_cross;
% eL_Sjj_cross = data_struct.eL_Sjj_cross;
% lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;
% XA_Sjj_cross = data_struct.XA_Sjj_cross;
% data_struct2.Svv_cross=Svv_cross;
% data_struct2.Sjj_cross=Sjj_cross;
% data_struct2.eL_Sjj_cross=eL_Sjj_cross;
% data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
% data_struct2.XA_Sjj_cross=XA_Sjj_cross;

% %%
% Define your parameters
Nr = size(Sjj_cross, 1); % Number of regions
Nw = Nw;%size(Sjj_cross, 3); % Number of frequencies
%Nsim = 50; % Number of simulations

% Precompute the upper triangular indices (excluding the diagonal)
idx_triu = triu(true(Nr), 1); % Only needs to be calculated once

% Preallocate memory for scores and distances based on estimated sizes
num_elements = sum(idx_triu(:)); % Number of upper triangular elements
total_size = Nsim * Nw * num_elements;

all_labels = zeros(total_size, 1);
all_eL_scores = zeros(total_size, 1);
all_lcmv_scores = zeros(total_size, 1);
all_XA_scores = zeros(total_size, 1);
%all_higgs_scores = zeros(total_size, 1);


% Preallocate distances
eL_distances_fro = zeros(Nsim * Nw, 1);
lcmv_distances_fro = zeros(Nsim * Nw, 1);
XA_distances_fro = zeros(Nsim * Nw, 1);
%higgs_distances_fro = zeros(Nsim * Nw, 1);


eL_distances_relative = zeros(Nsim * Nw, 1);
lcmv_distances_relative = zeros(Nsim * Nw, 1);
XA_distances_relative = zeros(Nsim * Nw, 1);
%higgs_distances_relative = zeros(Nsim * Nw, 1);


eL_distances_corr = zeros(Nsim * Nw, 1);
lcmv_distances_corr = zeros(Nsim * Nw, 1);
XA_distances_corr = zeros(Nsim * Nw, 1);
%higgs_distances_corr = zeros(Nsim * Nw, 1);

eL_distances_l1 = zeros(Nsim * Nw, 1); 
lcmv_distances_l1 = zeros(Nsim * Nw, 1);
XA_distances_l1 = zeros(Nsim * Nw, 1);
%higgs_distances_l1 = zeros(Nsim * Nw, 1);

% Counters for storing data efficiently
score_counter = 1;
distance_counter = 1;

for n = 1:10
    disp(n); % Display the simulation number
    for f = 1:Nw
        % Extract true and estimated covariance matrices (pre-compute abs)
        Sjj_true = abs((Sjj_cross(:, :, f, n)));
        eL_Sjj_est = abs((eL_Sjj_cross(:, :, f, n)));
        lcmv_Sjj_est = abs((lcmv_Sjj_cross(:, :, f, n)));
        XA_Sjj_est = abs((XA_Sjj_cross(:, :, f, n)));
       % higgs_Sjj_est = angle((higgs_Sjj_cross(:, :, f, n)));

        % Flatten the upper triangular parts of the matrices
        Sjj_true_vec = (Sjj_true(idx_triu));
        eL_Sjj_est_vec = (eL_Sjj_est(idx_triu));
        lcmv_Sjj_est_vec = (lcmv_Sjj_est(idx_triu));
        XA_Sjj_est_vec = (XA_Sjj_est(idx_triu));
        %higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));

        % Threshold for ground truth labels
        threshold = mean(Sjj_true_vec);
        labels = Sjj_true_vec > threshold;

        % Store labels and scores
        num_elements = numel(Sjj_true_vec);
        all_labels(score_counter:score_counter + num_elements - 1) = labels;
        all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
        all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
        all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
       % all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
        score_counter = score_counter + num_elements;

        % Compute Frobenius norm and relative error
        eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
        lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
        XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
        %higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');

        eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        %higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');

        % Compute correlation distances
        eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
        lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
        XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
        %higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
        % L1 norm
        eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
        lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
        XA_distances_l1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
        %higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
        distance_counter = distance_counter + 1;
    end
end


% Compute ROC curves using all data
[eL_X, eL_Y, ~, eL_AUC] = perfcurve(all_labels, all_eL_scores, true);
[lcmv_X, lcmv_Y, ~, lcmv_AUC] = perfcurve(all_labels, all_lcmv_scores, true);
[XA_X, XA_Y, ~, XA_AUC] = perfcurve(all_labels, all_XA_scores, true);
%[higgs_X, higgs_Y, ~, higgs_AUC] = perfcurve(all_labels, all_higgs_scores, true);




%here
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*
% Prepare data for daviolinplot of distances
methods = {'eLORETA', 'LCMV', 'XA'};%, 'higgs'};

% Ensure data are column vectors and apply the log(1 + x) transformation
fro_distances = {log(1 + eL_distances_fro(:)), log(1 + lcmv_distances_fro(:)), log(1 + XA_distances_fro(:))};%, log(1 + higgs_distances_fro(:))};
relative_distances = {log(1 + eL_distances_relative(:)), log(1 + lcmv_distances_relative(:)), log(1 + XA_distances_relative(:))};%, log(1 + higgs_distances_relative(:))};
corr_distances = {log(1+ eL_distances_corr(:)),log(1+lcmv_distances_corr(:)), log(1+XA_distances_corr(:))};%, log(1+higgs_distances_corr(:))};
l1_distances = {log(1+eL_distances_l1(:)), log(1+lcmv_distances_l1(:)),log(1+ XA_distances_l1(:))};%, log(1+ higgs_distances_l1(:))};

% Define colors for the methods (adjust these RGB values as desired)
colors = [0.2 0.6 0.8;   % Color for eLORETA
          0.8 0.4 0.2;   % Color for LCMV
          0.6 0.8 0.2;    % Color for XA
          0.5 0.5 0.5];  %Color higgs

% Plot daviolinplot for Frobenius norm distances without scatter data
% Prepare a new figure for all subplots in a 1x4 layout (excluding Correlation Distance)
figure;

% Define font sizes
titleFontSize = 18;
labelFontSize = 14;
legendFontSize = 8;

% Subplot 1: ROC Curve
subplot(1, 5, 1);
plot(eL_X, eL_Y, 'Color', colors(1,:), 'LineWidth', 2); hold on; % eLORETA
plot(lcmv_X, lcmv_Y, 'Color', colors(2,:), 'LineWidth', 2);      % LCMV
plot(XA_X, XA_Y, 'Color', colors(3,:), 'LineWidth', 2);          % XA
%plot(higgs_X, higgs_Y, 'Color', colors(4,:), 'LineWidth', 2);    % Higgs
legend(['eLORETA (AUC = ' num2str(eL_AUC, '%.2f') ')'], ...
      ['LCMV (AUC = ' num2str(lcmv_AUC, '%.2f') ')'], ...
      ['XA (AUC = ' num2str(XA_AUC, '%.2f') ')'], 'Location', 'Best', ...
      'FontSize', 10, 'FontWeight', 'bold');  % Set legend to bold and size 10
xlabel('False Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');  % Set axes labels to bold and size 14
ylabel('True Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');
title('Combined ROC', 'FontSize', 18, 'FontWeight', 'bold');  % Title to bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks and y-ticks to bold size 14
hold off;

% Subplot 2: Frobenius Norm Distance
subplot(1, 5, 2);
daviolinplot(fro_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('Frobenius Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('Frobenius Norm', 'FontSize', 18, 'FontWeight', 'bold');       % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 3: Relative Frobenius Error
subplot(1, 5, 3);
daviolinplot(relative_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('Frob Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('Frob Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 4: L1 Distance
subplot(1, 5, 4);
daviolinplot(l1_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('L1 Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('L1 Norm', 'FontSize', 18, 'FontWeight', 'bold');        % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 5: Spectral Norm Distance (L2 Distance)
subplot(1, 5, 5);
daviolinplot(corr_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('L1 Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('L1 Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Adjust the layout for better visualization
set(gcf, 'Position', [100, 100, 1800, 500]); % Adjust figure size for the 1x5 wide layout
set(gcf,'Color','w');

%% Read Data 
% Svv_cross = data_struct.Svv_cross;
% Sjj_cross = data_struct.Sjj_cross;
% eL_Sjj_cross = data_struct.eL_Sjj_cross;
% lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;
% XA_Sjj_cross = data_struct.XA_Sjj_cross;
% data_struct2.Svv_cross=Svv_cross;
% data_struct2.Sjj_cross=Sjj_cross;
% data_struct2.eL_Sjj_cross=eL_Sjj_cross;
% data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
% data_struct2.XA_Sjj_cross=XA_Sjj_cross;

% %%
% Define your parameters
Nr = size(Sjj_cross, 1); % Number of regions
Nw = Nw;%size(Sjj_cross, 3); % Number of frequencies
%Nsim = 50; % Number of simulations

% Precompute the upper triangular indices (excluding the diagonal)
idx_triu = triu(true(Nr), 1); % Only needs to be calculated once

% Preallocate memory for scores and distances based on estimated sizes
num_elements = sum(idx_triu(:)); % Number of upper triangular elements
total_size = Nsim * Nw * num_elements;

all_labels = zeros(total_size, 1);
all_eL_scores = zeros(total_size, 1);
all_lcmv_scores = zeros(total_size, 1);
all_XA_scores = zeros(total_size, 1);
%all_higgs_scores = zeros(total_size, 1);


% Preallocate distances
eL_distances_fro = zeros(Nsim * Nw, 1);
lcmv_distances_fro = zeros(Nsim * Nw, 1);
XA_distances_fro = zeros(Nsim * Nw, 1);
%higgs_distances_fro = zeros(Nsim * Nw, 1);


eL_distances_relative = zeros(Nsim * Nw, 1);
lcmv_distances_relative = zeros(Nsim * Nw, 1);
XA_distances_relative = zeros(Nsim * Nw, 1);
%higgs_distances_relative = zeros(Nsim * Nw, 1);


eL_distances_corr = zeros(Nsim * Nw, 1);
lcmv_distances_corr = zeros(Nsim * Nw, 1);
XA_distances_corr = zeros(Nsim * Nw, 1);
%higgs_distances_corr = zeros(Nsim * Nw, 1);

eL_distances_l1 = zeros(Nsim * Nw, 1); 
lcmv_distances_l1 = zeros(Nsim * Nw, 1);
XA_distances_l1 = zeros(Nsim * Nw, 1);
%higgs_distances_l1 = zeros(Nsim * Nw, 1);

% Counters for storing data efficiently
score_counter = 1;
distance_counter = 1;

for n = 1:10
    disp(n); % Display the simulation number
    for f = 1:Nw
        % Extract true and estimated covariance matrices (pre-compute abs)
        Sjj_true = angle((Sjj_cross(:, :, f, n)));
        eL_Sjj_est = angle((eL_Sjj_cross(:, :, f, n)));
        lcmv_Sjj_est = angle((lcmv_Sjj_cross(:, :, f, n)));
        XA_Sjj_est = angle((XA_Sjj_cross(:, :, f, n)));
       % higgs_Sjj_est = angle((higgs_Sjj_cross(:, :, f, n)));

        % Flatten the upper triangular parts of the matrices
        Sjj_true_vec = (Sjj_true(idx_triu));
        eL_Sjj_est_vec = (eL_Sjj_est(idx_triu));
        lcmv_Sjj_est_vec = (lcmv_Sjj_est(idx_triu));
        XA_Sjj_est_vec = (XA_Sjj_est(idx_triu));
        %higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));

        % Threshold for ground truth labels
        threshold = mean(Sjj_true_vec);
        labels = Sjj_true_vec > threshold;

        % Store labels and scores
        num_elements = numel(Sjj_true_vec);
        all_labels(score_counter:score_counter + num_elements - 1) = labels;
        all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
        all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
        all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
       % all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
        score_counter = score_counter + num_elements;

        % Compute Frobenius norm and relative error
        eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
        lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
        XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
        %higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');

        eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        %higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');

        % Compute correlation distances
        eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
        lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
        XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
        %higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
        % L1 norm
        eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
        lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
        XA_distances_l1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
        %higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
        distance_counter = distance_counter + 1;
    end
end


% Compute ROC curves using all data
[eL_X, eL_Y, ~, eL_AUC] = perfcurve(all_labels, all_eL_scores, true);
[lcmv_X, lcmv_Y, ~, lcmv_AUC] = perfcurve(all_labels, all_lcmv_scores, true);
[XA_X, XA_Y, ~, XA_AUC] = perfcurve(all_labels, all_XA_scores, true);
%[higgs_X, higgs_Y, ~, higgs_AUC] = perfcurve(all_labels, all_higgs_scores, true);




%here
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*
% Prepare data for daviolinplot of distances
methods = {'eLORETA', 'LCMV', 'XA'};%, 'higgs'};

% Ensure data are column vectors and apply the log(1 + x) transformation
fro_distances = {log(1 + eL_distances_fro(:)), log(1 + lcmv_distances_fro(:)), log(1 + XA_distances_fro(:))};%, log(1 + higgs_distances_fro(:))};
relative_distances = {log(1 + eL_distances_relative(:)), log(1 + lcmv_distances_relative(:)), log(1 + XA_distances_relative(:))};%, log(1 + higgs_distances_relative(:))};
corr_distances = {log(1+ eL_distances_corr(:)),log(1+lcmv_distances_corr(:)), log(1+XA_distances_corr(:))};%, log(1+higgs_distances_corr(:))};
l1_distances = {log(1+eL_distances_l1(:)), log(1+lcmv_distances_l1(:)),log(1+ XA_distances_l1(:))};%, log(1+ higgs_distances_l1(:))};

% Define colors for the methods (adjust these RGB values as desired)
colors = [0.2 0.6 0.8;   % Color for eLORETA
          0.8 0.4 0.2;   % Color for LCMV
          0.6 0.8 0.2;    % Color for XA
          0.5 0.5 0.5];  %Color higgs

% Plot daviolinplot for Frobenius norm distances without scatter data
% Prepare a new figure for all subplots in a 1x4 layout (excluding Correlation Distance)
figure;

% Define font sizes
titleFontSize = 18;
labelFontSize = 14;
legendFontSize = 8;

% Subplot 1: ROC Curve
subplot(1, 5, 1);
plot(eL_X, eL_Y, 'Color', colors(1,:), 'LineWidth', 2); hold on; % eLORETA
plot(lcmv_X, lcmv_Y, 'Color', colors(2,:), 'LineWidth', 2);      % LCMV
plot(XA_X, XA_Y, 'Color', colors(3,:), 'LineWidth', 2);          % XA
%plot(higgs_X, higgs_Y, 'Color', colors(4,:), 'LineWidth', 2);    % Higgs
legend(['eLORETA (AUC = ' num2str(eL_AUC, '%.2f') ')'], ...
      ['LCMV (AUC = ' num2str(lcmv_AUC, '%.2f') ')'], ...
      ['XA (AUC = ' num2str(XA_AUC, '%.2f') ')'], 'Location', 'Best', ...
      'FontSize', 10, 'FontWeight', 'bold');  % Set legend to bold and size 10
xlabel('False Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');  % Set axes labels to bold and size 14
ylabel('True Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');
title('Combined ROC', 'FontSize', 18, 'FontWeight', 'bold');  % Title to bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks and y-ticks to bold size 14
hold off;

% Subplot 2: Frobenius Norm Distance
subplot(1, 5, 2);
daviolinplot(fro_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('Frobenius Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('Frobenius Norm', 'FontSize', 18, 'FontWeight', 'bold');       % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 3: Relative Frobenius Error
subplot(1, 5, 3);
daviolinplot(relative_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('Frob Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('Frob Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 4: L1 Distance
subplot(1, 5, 4);
daviolinplot(l1_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('L1 Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('L1 Norm', 'FontSize', 18, 'FontWeight', 'bold');        % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Subplot 5: Spectral Norm Distance (L2 Distance)
subplot(1, 5, 5);
daviolinplot(corr_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
             'xtlabels', methods, 'scatter', 0, ...
             'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
             'boxspacing', 1.2); % Increased boxplot spacing
ylabel('L1 Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
title('L1 Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14

% Adjust the layout for better visualization
set(gcf, 'Position', [100, 100, 1800, 500]); % Adjust figure size for the 1x5 wide layout
set(gcf,'Color','w');

% 
% %%
% % Display mean and standard deviation of distances
% mean_eL_fro = mean(eL_distances_fro);
% std_eL_fro = std(eL_distances_fro);
% 
% mean_lcmv_fro = mean(lcmv_distances_fro);
% std_lcmv_fro = std(lcmv_distances_fro);
% 
% mean_XA_fro = mean(XA_distances_fro);
% std_XA_fro = std(XA_distances_fro);
% 
% fprintf('Frobenius Norm Distances:\n');
% fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_fro, std_eL_fro);
% fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_fro, std_lcmv_fro);
% fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_fro, std_XA_fro);
% 
% % Similarly for spectral norm
% mean_eL_spectral = mean(eL_distances_spectral);
% std_eL_spectral = std(eL_distances_spectral);
% 
% mean_lcmv_spectral = mean(lcmv_distances_spectral);
% std_lcmv_spectral = std(lcmv_distances_spectral);
% 
% mean_XA_spectral = mean(XA_distances_spectral);
% std_XA_spectral = std(XA_distances_spectral);
% 
% fprintf('\nSpectral Norm Distances:\n');
% fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_spectral, std_eL_spectral);
% fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_spectral, std_lcmv_spectral);
% fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_spectral, std_XA_spectral);
% 
% % And for relative error
% mean_eL_relative = mean(eL_distances_relative);
% std_eL_relative = std(eL_distances_relative);
% 
% mean_lcmv_relative = mean(lcmv_distances_relative);
% std_lcmv_relative = std(lcmv_distances_relative);
% 
% mean_XA_relative = mean(XA_distances_relative);
% std_XA_relative = std(XA_distances_relative);
% 
% fprintf('\nRelative Errors:\n');
% fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_relative, std_eL_relative);
% fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_relative, std_lcmv_relative);
% fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_relative, std_XA_relative);
% 
% % Display for correlation distance
% mean_eL_corr = mean(eL_distances_corr);
% std_eL_corr = std(eL_distances_corr);
% 
% mean_lcmv_corr = mean(lcmv_distances_corr);
% std_lcmv_corr = std(lcmv_distances_corr);
% 
% mean_XA_corr = mean(XA_distances_corr);
% std_XA_corr = std(XA_distances_corr);
% 
% fprintf('\nCorrelation Distances:\n');
% fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_corr, std_eL_corr);
% fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_corr, std_lcmv_corr);
% fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_corr, std_XA_corr);


%
% R = voxel_roi_map;
% preprocessing_velocity;
% K= parameters.Model.K;
% L=K;
% T = parameters.Model.T;
% K =pinv(K);
% U_map =K;
% for j=1:Nw
%      T_omega(:,:,j) = U_map*T(:,:,j);
% end
% parameters.Model.U= T_omega;
% load('Data/Average_Velocity_ROISpace/GPfit_Delay_Mean.mat');
% %
% x_mean = mean(x_avg,2);
% x_std  = std(x_avg',1)';
% [e_mean,a_mean,s2_mean]=x2v(x_mean);
% [e_std,a_std,s2_std]=x2v(x_std);
% 
% e_mean =max(e_mean).*ones(size(e_mean));
% a_mean =max(a_mean).*ones(size(a_mean));
% s2_mean =max(s2_mean).*ones(size(s2_mean));
% x_mean = v2x(e_mean,a_mean,s2_mean);
% 
% e_std =max(e_std).*ones(size(e_std));
% a_std =max(a_std).*ones(size(a_std));
% s2_std =max(s2_std).*ones(size(s2_std));
% x_std = v2x(e_std,a_mean,s2_std);
%%
 S = data_struct.CrossM;
 parameters.Model = Compact_Model;
parameters.Compact_Model = Compact_Model;
parameters.Dimensions = Dimensions;
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
parameters.Dimensions.Nv = Nr;
Nv = parameters.Dimensions.Nv; % Set voxel to roi
Nw = parameters.Dimensions.Nw;

import functions.*
import functions.auxx.Simulations.*
import functions.auxx.Simulations.mn_cross.*
import functions.auxx.DataPreprosessing.*
S = aveReference(S);
% S = J;
%S = Source_PSD;
freq = data_struct.freqrange(1:47);
l = zeros(size(S,1),47);
 Sjj = mn_cross(S, parameters);
for  i = 1:19
    for j = 1:47
        l(i,j) = (diag((Sjj(i,i,j))));
    end
end
plot(freq(1:47),mean(l',2))
%%
S = data_struct.CrossM;
freq = data_struct.freqrange;
l = zeros(size(S,1), 90);

% Calculate the l values
for i = 1:19
    for j = 1:90
        l(i,j) = 10 * log10(diag(S(i,i,j)));
    end
end

% Increase smoothness by increasing the number of interpolation points
x = 1:90; % Original x-values (data points corresponding to the frequency range)
xq = linspace(1, 90, 500); % Increased number of points for smoother curve

% Interpolate each row (i.e., each 'i' value) using spline interpolation
l_smooth = zeros(size(l, 1), 500); % Pre-allocate space for smoothed data
for i = 1:19
    l_smooth(i, :) = spline(x, l(i, :), xq); % Interpolate each row
end

% Interpolate the frequency range
freq_smooth = interp1(x, freq(1:90), xq, 'spline'); % Interpolate the frequency range using spline interpolation

% Plot the smoothed data with increased line width
figure; % Create a new figure window
plot(freq_smooth, l_smooth', 'LineWidth', 2) % Set LineWidth to 2 for thicker lines
xlabel('Frequency (Hz)');
ylabel('Log Power (dB)');
title('Smoothed Logarithmic Values of S');
%%
S = l;
[N_c, Nw] = size(l);

% Define Xi-Omega and Alpha-Omega model functions
xi_omega_model = @(e1, e2, e3, omega) e1 ./ (1 + e2 * omega.^2).^e3;
alpha_omega_model = @(a1, a2, a3, a4, omega) a1 ./ (1 + a2 * (omega - a3).^2).^a4;

% Set up the frequency vector (omega)
omega = freq; % Adjust the frequency range and number of points (Nw)

% Initialize storage for fitted parameters
xi_params = zeros(N_c, 3); % xi_omega parameters for each channel (e1, e2, e3)
alpha_params = zeros(N_c, 4); % alpha_omega parameters for each channel (a1, a2, a3, a4)

% Loop over each channel and perform fitting
for i = 1:N_c
    % Extract the spectrum for the current channel
    spectrum = S(i, :); % Spectrum for channel i, of size Nw
    
    % Define the combined model function for fitting (7 parameters total)
    combined_model = @(params, omega) ...
        xi_omega_model(params(1), params(2), params(3), omega) + ...
        alpha_omega_model(params(4), params(5), params(6), params(7), omega);
    
    % Define the error function (sum of squared errors)
    % Add penalty for a3 outside the interval [8, 13]
    error_func = @(params) sum((combined_model(params, omega) - spectrum).^2) + ...
        1e6 * (params(6) < 8 || params(6) > 13); % Penalize if a3 (params(6)) is outside [8, 13]
    
    % Set initial guess for the parameters (7 parameters total)
    initial_params = [1, 1, 1, 1, 1, 10, 1]; % Example initial guess, set a3 within the range [8, 13]
    
    % Transform parameters to ensure positivity (exponential transformation)
    transform_func = @(x) exp(x); % Exponentiate to ensure positivity
    
    % Use fminsearch to minimize the error function
    options = optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000);
    
    % Minimize the error function with transformed parameters
    fminsearch_error_func = @(transformed_params) error_func(transform_func(transformed_params));
    
    fitted_params_transformed = fminsearch(fminsearch_error_func, log(initial_params), options); % Log to apply exp later
    
    % Transform the fitted parameters back to positive values
    fitted_params = transform_func(fitted_params_transformed);
    
    % Store the fitted parameters for Xi and Alpha
    xi_params(i, :) = fitted_params(1:3); % Store Xi-Omega parameters (e1, e2, e3)
    alpha_params(i, :) = fitted_params(4:7); % Store Alpha-Omega parameters (a1, a2, a3, a4)
    
    % Plot the original spectrum and the fitted model
    figure;
    plot(omega, spectrum, 'b', 'LineWidth', 1.5); % Original spectrum in blue
    hold on;
    plot(omega, combined_model(fitted_params, omega), 'r--', 'LineWidth', 1.5); % Fitted model in red dashed line
    title(['Channel ' num2str(i) ' Fitting']);
    xlabel('Frequency (omega)');
    ylabel('Spectrum Amplitude');
    legend('Original Spectrum', 'Fitted Model');
    grid on;
    hold off;
end

% xi_params and alpha_params will now hold the fitted parameters for each channel



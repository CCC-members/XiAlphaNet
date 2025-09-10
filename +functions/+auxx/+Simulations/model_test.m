% Povilas Karvelis (2024). daviolinplot - beautiful violin and raincloud plots
% (https://github.com/povilaskarvelis/DataViz/releases/tag/v3.2.4), GitHub. Retrieved September 20, 2024.

% Clear workspace and initialize
clear all;
clc;

% Load structural data
load("/mnt/Develop/Ronaldo/dev/Data/NewFolder/structural/parameters.mat");

% Set parameters for simulation in ROI space
parameters.Model = Compact_Model;
parameters.Compact_Model = Compact_Model;
parameters.Dimensions = Dimensions;
parameters.Dimensions.Nv = parameters.Dimensions.Nr;
% Clear temporary variables to save memory
clear Compact_Model Model Dimensions;

% Import necessary functions
import functions.*;
import functions.auxx.*;
import functions.auxx.BayesOptimization.*;
import functions.auxx.CrossSpectrum.*;
import functions.auxx.ModelVectorization.*;
import functions.auxx.Simulations.*;
import functions.auxx.TOperator.*;
import tools.*;
import functions.auxx.DataPreprocessing.*;


% Set up simulation parameters
Ne = parameters.Dimensions.Ne;  % Number of electrodes
Nr = parameters.Dimensions.Nr;  % Number of ROIs
Nv = parameters.Dimensions.Nr;  % Set voxel to ROI  
Nw = parameters.Dimensions.Nw;  % Number of frequency bins
Nsim = 3;  % Number of simulations
N_wishart = 1000;
conn_spec = norm(parameters.Model.C,'fro');

% Xi-AlphaNET properties
disp("--> Estimating source cross-spectrum");

% Set model parameters for Xi-AlphaNET estimation
properties.model_params.nFreqs = Nw;
properties.model_params.BayesIter_Delay = 30;
properties.model_params.BayesIter_Reg1 = 20;
properties.model_params.BayesIter_Reg2 = 100;
properties.model_params.Nrand1 = 10;
properties.model_params.Nrand2 = 50;
properties.model_params.delay.lambda_space_cd = [[0.4, 1.6]; [10^(-10), 1/conn_spec]];
properties.general_params.parallel.conn_delay = 1;
properties.model_params.stoch1 = 1;
properties.model_params.stoch2 = 1;
properties.model_params.tensor_field.default = 0;

% Import necessary functions for Xi-AlphaNET estimation
import app.*;
import app.functions.*;
import functions.*;
import functions.StochasticFISTA.*;
import functions.auxx.*;
import functions.auxx.BayesOptimization.*;
import functions.auxx.CrossSpectrum.*;
import functions.auxx.ModelVectorization.*;
import functions.auxx.Regularization.*;
import functions.auxx.TOperator.*;
import tools.*;
import functions.auxx.Simulations.*;
import functions.auxx.Simulations.inverse.*;
import functions.auxx.Simulations.private.*;

% Load model and transformation matrices
L = parameters.Model.K;  % Transformation matrix for cross-spectrum

% Set data directory for simulation
dir_data = '/mnt/Develop/Ronaldo/dev/Data/norms';
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
selected_folders = subject_folders(randperm(length(subject_folders), Nsim));

K = parameters.Model.K;
R = parameters.Compact_Model.R;
% Simulation loop

for j = 1:Nsim

    tic;
    disp(['Processing simulation ', num2str(j), ' of ', num2str(Nsim)]);
%%
    % Select subject folder and load data
    subject_folder = selected_folders(j).name;
    mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
    data_struct = load(mat_file_path);
    
    % Extract cross-spectrum and reference it
    disp('->> Reading Scalp Cross')
    Svv = data_struct.data_struct.CrossM(:,:,1:Nw);
    Svv = functions.auxx.DataPreprosessing.aveReference(Svv);

     % Set frequency range and parameters for the current simulation
    freq = data_struct.data_struct.freqrange(1:Nw);
    parameters.Data.freq = freq;
    %%
    % Source cross-spectrum 
    disp('->> Estimating Source Cross with MN')
    Sjj = mn_cross(Svv,K,0);

    % Generate and process simulated cross-spectrum
    disp('->> Simulating Scalp Cross Wishart Noise + FModel')
    parfor i = 1:Nw
        i
        % Generate complex Wishart matrices and simulate cross-spectrum
        Sjj_cross(:,:,i) = generate_complex_wishart(Sjj(:,:,i), N_wishart);
        Svv_cross(:,:,i) = L * Sjj_cross(:,:,i) * L';  % Apply transformation
    end
    
    toc;  % End of simulation processing
    
    
    % Prepare data for Xi-AlphaNET estimation
    data.Cross = Svv_cross;
    data.age = 25;  % Age of the subject
    data.freq = freq;
   
    % Perform Xi-AlphaNET estimation
    disp('->> Xi-AlphaNeT Inverse Solution')
    [x, ~, G,x0] = Xi_ALphaNET(properties, data, parameters);
    [source_act_cross] = functions.auxx.CrossSpectrum.eval_source_conn(x.Solution, data.freq, parameters.Model.R, properties, parameters);
    % Store the results of the Xi-AlphaNET simulation
    XA_Sjj_cross = source_act_cross.Cross.Full;
    
    % Parallelized loop for source estimation using different methods
    disp('->> eLORETA and LCMV Inverse Solutions')
    parfor i = 1:Nw
        % Perform eLORETA and LCMV source estimations
        source = inverse(Svv_cross(:,:,i), L);
        eL_Sjj_cross(:,:,i) = source.eloreata.Sjj;
        lcmv_Sjj_cross(:,:,i) = source.lcmv.Sjj;
    end
    toc;
  
    % Precompute coherence and phase values
    disp('->> Evaluating Distances')
    import functions.auxx.OptimizedOperations.*
    coherence_XA = coherence(XA_Sjj_cross);
    coherence_Sjj = coherence(Sjj_cross);
    angle_XA = angle(XA_Sjj_cross);
    angle_Sjj = angle(Sjj_cross);
    
    % Log-spectrum benchmarking for XA, eL, and LCMV
   
    for i = 1:Nw
        xl(:,i) = log(real(diag(XA_Sjj_cross(:,:,i))));
        el(:,i) = log(real(diag(eL_Sjj_cross(:,:,i))));
        cl(:,i) = log(real(diag(lcmv_Sjj_cross(:,:,i))));
        l(:,i) = log(real(diag(Sjj_cross(:,:,i))));
    end
   
    
    % Benchmark log-spectrum for XA, eL, and LCMV
    XA_Spectrum_Fro(j) = tensor_norm(xl - l, 2);
    XA_Spectrum_FroR(j) = XA_Spectrum_Fro(j) / tensor_norm(l, 2);
    XA_Spectrum_L1(j) = tensor_norm(xl - l, 1);
    XA_Spectrum_L1R(j) = XA_Spectrum_L1(j) / tensor_norm(l, 1);
    
    eL_Spectrum_Fro(j) = tensor_norm(el - l, 2);
    eL_Spectrum_FroR(j) = eL_Spectrum_Fro(j) / tensor_norm(l, 2);
    eL_Spectrum_L1(j) = tensor_norm(el - l, 1);
    eL_Spectrum_L1R(j) = eL_Spectrum_L1(j) / tensor_norm(l, 1);
    
    cl_Spectrum_Fro(j) = tensor_norm(cl - l, 2);
    cl_Spectrum_FroR(j) = cl_Spectrum_Fro(j) / tensor_norm(l, 2);
    cL_Spectrum_L1(j) = tensor_norm(cl - l, 1);
    cL_Spectrum_L1R(j) = cL_Spectrum_L1(j) / tensor_norm(l, 1);
    
    % Benchmark coherence XA, eL, and LCMV
    XA_Coherence_Fro(j) = tensor_norm(coherence_XA - coherence_Sjj, 2);
    XA_Coherence_FroR(j) = XA_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    XA_Coherence_L1(j) = tensor_norm(coherence_XA - coherence_Sjj, 1);
    XA_Coherence_L1R(j) = XA_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    
    eL_Coherence_Fro(j) = tensor_norm(coherence(eL_Sjj_cross) - coherence_Sjj, 2);
    eL_Coherence_FroR(j) = eL_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    eL_Coherence_L1(j) = tensor_norm(coherence(eL_Sjj_cross) - coherence_Sjj, 1);
    eL_Coherence_L1R(j) = eL_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    
    lcmv_Coherence_Fro(j) = tensor_norm(coherence(lcmv_Sjj_cross) - coherence_Sjj, 2);
    lcmv_Coherence_FroR(j) = lcmv_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    lcmv_Coherence_L1(j) = tensor_norm(coherence(lcmv_Sjj_cross) - coherence_Sjj, 1);
    lcmv_Coherence_L1R(j) = lcmv_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    
    % Benchmark phase XA, eL, and LCMV
    XA_Phase_Fro(j) = tensor_norm(angle_XA - angle_Sjj, 2);
    XA_Phase_FroR(j) = XA_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    XA_Phase_L1(j) = tensor_norm(angle_XA - angle_Sjj, 1);
    XA_Phase_L1R(j) = XA_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    
    eL_Phase_Fro(j) = tensor_norm(angle(eL_Sjj_cross) - angle_Sjj, 2);
    eL_Phase_FroR(j) = eL_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    eL_Phase_L1(j) = tensor_norm(angle(eL_Sjj_cross) - angle_Sjj, 1);
    eL_Phase_L1R(j) = eL_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    
    lcmv_Phase_Fro(j) = tensor_norm(angle(lcmv_Sjj_cross) - angle_Sjj, 2);
    lcmv_Phase_FroR(j) = lcmv_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    lcmv_Phase_L1(j) = tensor_norm(angle(lcmv_Sjj_cross) - angle_Sjj, 1);
    lcmv_Phase_L1R(j) = lcmv_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
   
end

%%
figure(1)
hold on
plot(freq,mean(xl',2),'Color','r','LineWidth',3)
plot(freq,mean(l',2),'-','Color','b')
hold off
norm(xl-l,'fro')/norm(l,'fro')
norm(el-l,'fro')/norm(l,'fro')


figure(2)
hold on
plot(freq,mean(el',2),'Color','r')
plot(freq,mean(l',2),'-','Color','b')
hold off
%%

%% === Log-Spectrum ===
% Relative Frobenius Norm
data_logspec_FroR = [XA_Spectrum_FroR(:); eL_Spectrum_FroR(:); cl_Spectrum_FroR(:)];
group_logspec_FroR = [ ...
    repmat({'XA'}, length(XA_Spectrum_FroR), 1);
    repmat({'eLORETA'}, length(eL_Spectrum_FroR), 1);
    repmat({'LCMV'}, length(cl_Spectrum_FroR), 1)];
figure;
boxplot(data_logspec_FroR, group_logspec_FroR);
ylabel('Relative Frobenius Norm');
title('Log-Spectrum Relative Frobenius Norm');
grid on;
hold on;
medians = [median(XA_Spectrum_FroR), median(eL_Spectrum_FroR), median(cl_Spectrum_FroR)];
[~, bestIdx] = min(medians);
ylims = ylim;
y_offset = 0.05 * (ylims(2) - ylims(1));
ylim([ylims(1), ylims(2) + y_offset]);
text(bestIdx, ylims(2) + 0.02 * (ylims(2) - ylims(1)), '*', 'FontSize', 20, 'HorizontalAlignment', 'center');

% Relative L1 Norm
data_logspec_L1R = [XA_Spectrum_L1R(:); eL_Spectrum_L1R(:); cL_Spectrum_L1R(:)];
group_logspec_L1R = [ ...
    repmat({'XA'}, length(XA_Spectrum_L1R), 1);
    repmat({'eLORETA'}, length(eL_Spectrum_L1R), 1);
    repmat({'LCMV'}, length(cL_Spectrum_L1R), 1)];
figure;
boxplot(data_logspec_L1R, group_logspec_L1R);
ylabel('Relative L1 Norm');
title('Log-Spectrum Relative L1 Norm');
grid on;
hold on;
medians = [median(XA_Spectrum_L1R), median(eL_Spectrum_L1R), median(cL_Spectrum_L1R)];
[~, bestIdx] = min(medians);
ylims = ylim;
y_offset = 0.05 * (ylims(2) - ylims(1));
ylim([ylims(1), ylims(2) + y_offset]);
text(bestIdx, ylims(2) + 0.02 * (ylims(2) - ylims(1)), '*', 'FontSize', 20, 'HorizontalAlignment', 'center');

%% === Coherence ===
% Relative Frobenius Norm
data_coh_FroR = [XA_Coherence_FroR(:); eL_Coherence_FroR(:); lcmv_Coherence_FroR(:)];
group_coh_FroR = [ ...
    repmat({'XA'}, length(XA_Coherence_FroR), 1);
    repmat({'eLORETA'}, length(eL_Coherence_FroR), 1);
    repmat({'LCMV'}, length(lcmv_Coherence_FroR), 1)];
figure;
boxplot(data_coh_FroR, group_coh_FroR);
ylabel('Relative Frobenius Norm');
title('Coherence Relative Frobenius Norm');
grid on;
hold on;
medians = [median(XA_Coherence_FroR), median(eL_Coherence_FroR), median(lcmv_Coherence_FroR)];
[~, bestIdx] = min(medians);
ylims = ylim;
y_offset = 0.05 * (ylims(2) - ylims(1));
ylim([ylims(1), ylims(2) + y_offset]);
text(bestIdx, ylims(2) + 0.02 * (ylims(2) - ylims(1)), '*', 'FontSize', 20, 'HorizontalAlignment', 'center');

% Relative L1 Norm
data_coh_L1R = [XA_Coherence_L1R(:); eL_Coherence_L1R(:); lcmv_Coherence_L1R(:)];
group_coh_L1R = [ ...
    repmat({'XA'}, length(XA_Coherence_L1R), 1);
    repmat({'eLORETA'}, length(eL_Coherence_L1R), 1);
    repmat({'LCMV'}, length(lcmv_Coherence_L1R), 1)];
figure;
boxplot(data_coh_L1R, group_coh_L1R);
ylabel('Relative L1 Norm');
title('Coherence Relative L1 Norm');
grid on;
hold on;
medians = [median(XA_Coherence_L1R), median(eL_Coherence_L1R), median(lcmv_Coherence_L1R)];
[~, bestIdx] = min(medians);
ylims = ylim;
y_offset = 0.05 * (ylims(2) - ylims(1));
ylim([ylims(1), ylims(2) + y_offset]);
text(bestIdx, ylims(2) + 0.02 * (ylims(2) - ylims(1)), '*', 'FontSize', 20, 'HorizontalAlignment', 'center');

%% === Phase ===
% Relative Frobenius Norm
data_phase_FroR = [XA_Phase_FroR(:); eL_Phase_FroR(:); lcmv_Phase_FroR(:)];
group_phase_FroR = [ ...
    repmat({'XA'}, length(XA_Phase_FroR), 1);
    repmat({'eLORETA'}, length(eL_Phase_FroR), 1);
    repmat({'LCMV'}, length(lcmv_Phase_FroR), 1)];
figure;
boxplot(data_phase_FroR, group_phase_FroR);
ylabel('Relative Frobenius Norm');
title('Phase Relative Frobenius Norm');
grid on;
hold on;
medians = [median(XA_Phase_FroR), median(eL_Phase_FroR), median(lcmv_Phase_FroR)];
[~, bestIdx] = min(medians);
ylims = ylim;
y_offset = 0.05 * (ylims(2) - ylims(1));
ylim([ylims(1), ylims(2) + y_offset]);
text(bestIdx, ylims(2) + 0.02 * (ylims(2) - ylims(1)), '*', 'FontSize', 20, 'HorizontalAlignment', 'center');

% Relative L1 Norm
data_phase_L1R = [XA_Phase_L1R(:); eL_Phase_L1R(:); lcmv_Phase_L1R(:)];
group_phase_L1R = [ ...
    repmat({'XA'}, length(XA_Phase_L1R), 1);
    repmat({'eLORETA'}, length(eL_Phase_L1R), 1);
    repmat({'LCMV'}, length(lcmv_Phase_L1R), 1)];
figure;
boxplot(data_phase_L1R, group_phase_L1R);
ylabel('Relative L1 Norm');
title('Phase Relative L1 Norm');
grid on;
hold on;
medians = [median(XA_Phase_L1R), median(eL_Phase_L1R), median(lcmv_Phase_L1R)];
[~, bestIdx] = min(medians);
ylims = ylim;
y_offset = 0.05 * (ylims(2) - ylims(1));
ylim([ylims(1), ylims(2) + y_offset]);
text(bestIdx, ylims(2) + 0.02 * (ylims(2) - ylims(1)), '*', 'FontSize', 20, 'HorizontalAlignment', 'center');
%%


% %% Read Data  and validate spectra
% % Svv_cross = data_struct.Svv_cross;
% % Sjj_cross = data_struct.Sjj_cross;
% % eL_Sjj_cross = data_struct.eL_Sjj_cross;
% % lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;c
% % XA_Sjj_cross = data_struct.XA_Sjj_cross;
% % data_struct2.Svv_cross=Svv_cross;
% % data_struct2.Sjj_cross=Sjj_cross;
% % data_struct2.eL_Sjj_cross=eL_Sjj_cross;
% % data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
% % data_struct2.XA_Sjj_cross=XA_Sjj_cross;
% 
% % %%
% % Define your parameters
% Nr = size(Sjj_cross, 1); % Number of regions
% Nw = Nw;%size(Sjj_cross, 3); % Number of frequencies
% %Nsim = 50; % Number of simulations
% 
% % Precompute the upper triangular indices (excluding the diagonal)
% idx_triu = triu(true(Nr), 1); % Only needs to be calculated once
% 
% % Preallocate memory for scores and distances based on estimated sizes
% num_elements = sum(idx_triu(:)); % Number of upper triangular elements
% total_size = Nsim * Nw * num_elements;
% 
% all_labels = zeros(total_size, 1);
% all_eL_scores = zeros(total_size, 1);
% all_lcmv_scores = zeros(total_size, 1);
% all_XA_scores = zeros(total_size, 1);
% %all_higgs_scores = zeros(total_size, 1);
% 
% 
% % Preallocate distances
% eL_distances_fro = zeros(Nsim * Nw, 1);
% lcmv_distances_fro = zeros(Nsim * Nw, 1);
% XA_distances_fro = zeros(Nsim * Nw, 1);
% %higgs_distances_fro = zeros(Nsim * Nw, 1);
% 
% 
% eL_distances_relative = zeros(Nsim * Nw, 1);
% lcmv_distances_relative = zeros(Nsim * Nw, 1);
% XA_distances_relative = zeros(Nsim * Nw, 1);
% %higgs_distances_relative = zeros(Nsim * Nw, 1);
% 
% 
% eL_distances_corr = zeros(Nsim * Nw, 1);
% lcmv_distances_corr = zeros(Nsim * Nw, 1);
% XA_distances_corr = zeros(Nsim * Nw, 1);
% %higgs_distances_corr = zeros(Nsim * Nw, 1);
% 
% eL_distances_l1 = zeros(Nsim * Nw, 1); 
% lcmv_distances_l1 = zeros(Nsim * Nw, 1);
% XA_distances_l1 = zeros(Nsim * Nw, 1);
% %higgs_distances_l1 = zeros(Nsim * Nw, 1);
% 
% % Counters for storing data efficiently
% score_counter = 1;
% distance_counter = 1;
% 
% for n = 1:1
%     disp(n); % Display the simulation number
%     for f = 1:Nw
%         % Extract true and estimated covariance matrices (pre-compute abs)
%         Sjj_true = (abs((Sjj_cross(:, :, f, n))));
%         eL_Sjj_est = (abs((eL_Sjj_cross(:, :, f, n))));
%         lcmv_Sjj_est = (abs((lcmv_Sjj_cross(:, :, f, n))));
%         XA_Sjj_est = (abs((XA_Sjj_cross(:, :, f, n))));
%        % higgs_Sjj_est = angle((higgs_Sjj_cross(:, :, f, n)));
% 
%         % Flatten the upper triangular parts of the matrices
%         Sjj_true_vec = diag(Sjj_true);%(idx_triu));
%         eL_Sjj_est_vec = diag(eL_Sjj_est);%(idx_triu));
%         lcmv_Sjj_est_vec = diag(lcmv_Sjj_est);%(idx_triu));
%         XA_Sjj_est_vec = diag(XA_Sjj_est);%(idx_triu));
%         %higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));
% 
%         % Threshold for ground truth labels
%         threshold = prctile(Sjj_true_vec,75);
%         labels = Sjj_true_vec > threshold;
% 
%         % Store labels and scores
%         num_elements = numel(Sjj_true_vec);
%         all_labels(score_counter:score_counter + num_elements - 1) = labels;
%         all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
%         all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
%         all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
%        % all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
%         score_counter = score_counter + num_elements;
% 
%         % Compute Frobenius norm and relative error
%         eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
%         lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
%         XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
%         %higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');
% 
%         eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         %higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
% 
%         % Compute correlation distances
%         eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
%         lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
%         XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
%         %higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
%         % L1 norm
%         eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
%         lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
%         XA_distances_l1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
%         %higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
%         distance_counter = distance_counter + 1;
%     end
% end
% 
% 
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
% 
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
% 
% % Subplot 3: Relative Frobenius Error
% subplot(1, 5, 3);
% daviolinplot(relative_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('Frob Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('Frob Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 4: L1 Distance
% subplot(1, 5, 4);
% daviolinplot(l1_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('L1 Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('L1 Norm', 'FontSize', 18, 'FontWeight', 'bold');        % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 5: Spectral Norm Distance (L2 Distance)
% subplot(1, 5, 5);
% daviolinplot(corr_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('L1 Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('L1 Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Adjust the layout for better visualization
% set(gcf, 'Position', [100, 100, 1800, 500]); % Adjust figure size for the 1x5 wide layout
% set(gcf,'Color','w');
% 
% %% Read Data 
% % Svv_cross = data_struct.Svv_cross;
% % Sjj_cross = data_struct.Sjj_cross;
% % eL_Sjj_cross = data_struct.eL_Sjj_cross;
% % lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;
% % XA_Sjj_cross = data_struct.XA_Sjj_cross;
% % data_struct2.Svv_cross=Svv_cross;
% % data_struct2.Sjj_cross=Sjj_cross;
% % data_struct2.eL_Sjj_cross=eL_Sjj_cross;
% % data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
% % data_struct2.XA_Sjj_cross=XA_Sjj_cross;
% 
% % %%
% % Define your parameters
% Nr = size(Sjj_cross, 1); % Number of regions
% Nw = Nw;%size(Sjj_cross, 3); % Number of frequencies
% %Nsim = 50; % Number of simulations
% 
% % Precompute the upper triangular indices (excluding the diagonal)
% idx_triu = triu(true(Nr), 1); % Only needs to be calculated once
% 
% % Preallocate memory for scores and distances based on estimated sizes
% num_elements = sum(idx_triu(:)); % Number of upper triangular elements
% total_size = Nsim * Nw * num_elements;
% 
% all_labels = zeros(total_size, 1);
% all_eL_scores = zeros(total_size, 1);
% all_lcmv_scores = zeros(total_size, 1);
% all_XA_scores = zeros(total_size, 1);
% %all_higgs_scores = zeros(total_size, 1);
% 
% 
% % Preallocate distances
% eL_distances_fro = zeros(Nsim * Nw, 1);
% lcmv_distances_fro = zeros(Nsim * Nw, 1);
% XA_distances_fro = zeros(Nsim * Nw, 1);
% %higgs_distances_fro = zeros(Nsim * Nw, 1);
% 
% 
% eL_distances_relative = zeros(Nsim * Nw, 1);
% lcmv_distances_relative = zeros(Nsim * Nw, 1);
% XA_distances_relative = zeros(Nsim * Nw, 1);
% %higgs_distances_relative = zeros(Nsim * Nw, 1);
% 
% 
% eL_distances_corr = zeros(Nsim * Nw, 1);
% lcmv_distances_corr = zeros(Nsim * Nw, 1);
% XA_distances_corr = zeros(Nsim * Nw, 1);
% %higgs_distances_corr = zeros(Nsim * Nw, 1);
% 
% eL_distances_l1 = zeros(Nsim * Nw, 1); 
% lcmv_distances_l1 = zeros(Nsim * Nw, 1);
% XA_distances_l1 = zeros(Nsim * Nw, 1);
% %higgs_distances_l1 = zeros(Nsim * Nw, 1);
% 
% % Counters for storing data efficiently
% score_counter = 1;
% distance_counter = 1;
% 
% for n = 1:1
%     disp(n); % Display the simulation number
%     for f = 1:Nw
%         % Extract true and estimated covariance matrices (pre-compute abs)
%         Sjj_true = abs((Sjj_cross(:, :, f, n)));
%         eL_Sjj_est = abs((eL_Sjj_cross(:, :, f, n)));
%         lcmv_Sjj_est = abs((lcmv_Sjj_cross(:, :, f, n)));
%         XA_Sjj_est = abs((XA_Sjj_cross(:, :, f, n)));
%        % higgs_Sjj_est = angle((higgs_Sjj_cross(:, :, f, n)));
% 
%         % Flatten the upper triangular parts of the matrices
%         Sjj_true_vec = (Sjj_true(idx_triu));
%         eL_Sjj_est_vec = (eL_Sjj_est(idx_triu));
%         lcmv_Sjj_est_vec = (lcmv_Sjj_est(idx_triu));
%         XA_Sjj_est_vec = (XA_Sjj_est(idx_triu));
%         %higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));
% 
%         % Threshold for ground truth labels
%         threshold = mean(Sjj_true_vec);
%         labels = Sjj_true_vec > threshold;
% 
%         % Store labels and scores
%         num_elements = numel(Sjj_true_vec);
%         all_labels(score_counter:score_counter + num_elements - 1) = labels;
%         all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
%         all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
%         all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
%        % all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
%         score_counter = score_counter + num_elements;
% 
%         % Compute Frobenius norm and relative error
%         eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
%         lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
%         XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
%         %higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');
% 
%         eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         %higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
% 
%         % Compute correlation distances
%         eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
%         lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
%         XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
%         %higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
%         % L1 norm
%         eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
%         lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
%         XA_distances_l1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
%         %higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
%         distance_counter = distance_counter + 1;
%     end
% end
% 
% 
% % Compute ROC curves using all data
% [eL_X, eL_Y, ~, eL_AUC] = perfcurve(all_labels, all_eL_scores, true);
% [lcmv_X, lcmv_Y, ~, lcmv_AUC] = perfcurve(all_labels, all_lcmv_scores, true);
% [XA_X, XA_Y, ~, XA_AUC] = perfcurve(all_labels, all_XA_scores, true);
% %[higgs_X, higgs_Y, ~, higgs_AUC] = perfcurve(all_labels, all_higgs_scores, true);
% 
% 
% 
% 
% %here
% import guide.Visualization.*
% import guide.Visualization.DataVizm.*
% import guide.Visualization.DataVizm.daviolinplot.*
% % Prepare data for daviolinplot of distances
% methods = {'eLORETA', 'LCMV', 'XA'};%, 'higgs'};
% 
% % Ensure data are column vectors and apply the log(1 + x) transformation
% fro_distances = {log(1 + eL_distances_fro(:)), log(1 + lcmv_distances_fro(:)), log(1 + XA_distances_fro(:))};%, log(1 + higgs_distances_fro(:))};
% relative_distances = {log(1 + eL_distances_relative(:)), log(1 + lcmv_distances_relative(:)), log(1 + XA_distances_relative(:))};%, log(1 + higgs_distances_relative(:))};
% corr_distances = {log(1+ eL_distances_corr(:)),log(1+lcmv_distances_corr(:)), log(1+XA_distances_corr(:))};%, log(1+higgs_distances_corr(:))};
% l1_distances = {log(1+eL_distances_l1(:)), log(1+lcmv_distances_l1(:)),log(1+ XA_distances_l1(:))};%, log(1+ higgs_distances_l1(:))};
% 
% % Define colors for the methods (adjust these RGB values as desired)
% colors = [0.2 0.6 0.8;   % Color for eLORETA
%           0.8 0.4 0.2;   % Color for LCMV
%           0.6 0.8 0.2;    % Color for XA
%           0.5 0.5 0.5];  %Color higgs
% 
% % Plot daviolinplot for Frobenius norm distances without scatter data
% % Prepare a new figure for all subplots in a 1x4 layout (excluding Correlation Distance)
% figure;
% 
% % Define font sizes
% titleFontSize = 18;
% labelFontSize = 14;
% legendFontSize = 8;
% 
% % Subplot 1: ROC Curve
% subplot(1, 5, 1);
% plot(eL_X, eL_Y, 'Color', colors(1,:), 'LineWidth', 2); hold on; % eLORETA
% plot(lcmv_X, lcmv_Y, 'Color', colors(2,:), 'LineWidth', 2);      % LCMV
% plot(XA_X, XA_Y, 'Color', colors(3,:), 'LineWidth', 2);          % XA
% %plot(higgs_X, higgs_Y, 'Color', colors(4,:), 'LineWidth', 2);    % Higgs
% legend(['eLORETA (AUC = ' num2str(eL_AUC, '%.2f') ')'], ...
%       ['LCMV (AUC = ' num2str(lcmv_AUC, '%.2f') ')'], ...
%       ['XA (AUC = ' num2str(XA_AUC, '%.2f') ')'], 'Location', 'Best', ...
%       'FontSize', 10, 'FontWeight', 'bold');  % Set legend to bold and size 10
% xlabel('False Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');  % Set axes labels to bold and size 14
% ylabel('True Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');
% title('Combined ROC', 'FontSize', 18, 'FontWeight', 'bold');  % Title to bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks and y-ticks to bold size 14
% hold off;
% 
% % Subplot 2: Frobenius Norm Distance
% subplot(1, 5, 2);
% daviolinplot(fro_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('Frobenius Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('Frobenius Norm', 'FontSize', 18, 'FontWeight', 'bold');       % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 3: Relative Frobenius Error
% subplot(1, 5, 3);
% daviolinplot(relative_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('Frob Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('Frob Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 4: L1 Distance
% subplot(1, 5, 4);
% daviolinplot(l1_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('L1 Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('L1 Norm', 'FontSize', 18, 'FontWeight', 'bold');        % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 5: Spectral Norm Distance (L2 Distance)
% subplot(1, 5, 5);
% daviolinplot(corr_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('L1 Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('L1 Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Adjust the layout for better visualization
% set(gcf, 'Position', [100, 100, 1800, 500]); % Adjust figure size for the 1x5 wide layout
% set(gcf,'Color','w');
% 
% %% Read Data 
% % Svv_cross = data_struct.Svv_cross;
% % Sjj_cross = data_struct.Sjj_cross;
% % eL_Sjj_cross = data_struct.eL_Sjj_cross;
% % lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;
% % XA_Sjj_cross = data_struct.XA_Sjj_cross;
% % data_struct2.Svv_cross=Svv_cross;
% % data_struct2.Sjj_cross=Sjj_cross;
% % data_struct2.eL_Sjj_cross=eL_Sjj_cross;
% % data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
% % data_struct2.XA_Sjj_cross=XA_Sjj_cross;
% 
% % %%
% % Define your parameters
% Nr = size(Sjj_cross, 1); % Number of regions
% Nw = Nw;%size(Sjj_cross, 3); % Number of frequencies
% %Nsim = 50; % Number of simulations
% 
% % Precompute the upper triangular indices (excluding the diagonal)
% idx_triu = triu(true(Nr), 1); % Only needs to be calculated once
% 
% % Preallocate memory for scores and distances based on estimated sizes
% num_elements = sum(idx_triu(:)); % Number of upper triangular elements
% total_size = Nsim * Nw * num_elements;
% 
% all_labels = zeros(total_size, 1);
% all_eL_scores = zeros(total_size, 1);
% all_lcmv_scores = zeros(total_size, 1);
% all_XA_scores = zeros(total_size, 1);
% %all_higgs_scores = zeros(total_size, 1);
% 
% 
% % Preallocate distances
% eL_distances_fro = zeros(Nsim * Nw, 1);
% lcmv_distances_fro = zeros(Nsim * Nw, 1);
% XA_distances_fro = zeros(Nsim * Nw, 1);
% %higgs_distances_fro = zeros(Nsim * Nw, 1);
% 
% 
% eL_distances_relative = zeros(Nsim * Nw, 1);
% lcmv_distances_relative = zeros(Nsim * Nw, 1);
% XA_distances_relative = zeros(Nsim * Nw, 1);
% %higgs_distances_relative = zeros(Nsim * Nw, 1);
% 
% 
% eL_distances_corr = zeros(Nsim * Nw, 1);
% lcmv_distances_corr = zeros(Nsim * Nw, 1);
% XA_distances_corr = zeros(Nsim * Nw, 1);
% %higgs_distances_corr = zeros(Nsim * Nw, 1);
% 
% eL_distances_l1 = zeros(Nsim * Nw, 1); 
% lcmv_distances_l1 = zeros(Nsim * Nw, 1);
% XA_distances_l1 = zeros(Nsim * Nw, 1);
% %higgs_distances_l1 = zeros(Nsim * Nw, 1);
% 
% % Counters for storing data efficiently
% score_counter = 1;
% distance_counter = 1;
% 
% for n = 1:1
%     disp(n); % Display the simulation number
%     for f = 1:Nw
%         % Extract true and estimated covariance matrices (pre-compute abs)
%         Sjj_true = angle((Sjj_cross(:, :, f, n)));
%         eL_Sjj_est = angle((eL_Sjj_cross(:, :, f, n)));
%         lcmv_Sjj_est = angle((lcmv_Sjj_cross(:, :, f, n)));
%         XA_Sjj_est = angle((XA_Sjj_cross(:, :, f, n)));
%        % higgs_Sjj_est = angle((higgs_Sjj_cross(:, :, f, n)));
% 
%         % Flatten the upper triangular parts of the matrices
%         Sjj_true_vec = (Sjj_true(idx_triu));
%         eL_Sjj_est_vec = (eL_Sjj_est(idx_triu));
%         lcmv_Sjj_est_vec = (lcmv_Sjj_est(idx_triu));
%         XA_Sjj_est_vec = (XA_Sjj_est(idx_triu));
%         %higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));
% 
%         % Threshold for ground truth labels
%         threshold = mean(Sjj_true_vec);
%         labels = Sjj_true_vec > threshold;
% 
%         % Store labels and scores
%         num_elements = numel(Sjj_true_vec);
%         all_labels(score_counter:score_counter + num_elements - 1) = labels;
%         all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
%         all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
%         all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
%        % all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
%         score_counter = score_counter + num_elements;
% 
%         % Compute Frobenius norm and relative error
%         eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
%         lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
%         XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
%         %higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');
% 
%         eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
%         %higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
% 
%         % Compute correlation distances
%         eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
%         lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
%         XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
%         %higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
%         % L1 norm
%         eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
%         lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
%         XA_distances_l1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
%         %higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
%         distance_counter = distance_counter + 1;
%     end
% end
% 
% 
% % Compute ROC curves using all data
% [eL_X, eL_Y, ~, eL_AUC] = perfcurve(all_labels, all_eL_scores, true);
% [lcmv_X, lcmv_Y, ~, lcmv_AUC] = perfcurve(all_labels, all_lcmv_scores, true);
% [XA_X, XA_Y, ~, XA_AUC] = perfcurve(all_labels, all_XA_scores, true);
% %[higgs_X, higgs_Y, ~, higgs_AUC] = perfcurve(all_labels, all_higgs_scores, true);
% 
% 
% 
% 
% %here
% import guide.Visualization.*
% import guide.Visualization.DataVizm.*
% import guide.Visualization.DataVizm.daviolinplot.*
% % Prepare data for daviolinplot of distances
% methods = {'eLORETA', 'LCMV', 'XA'};%, 'higgs'};
% 
% % Ensure data are column vectors and apply the log(1 + x) transformation
% fro_distances = {log(1 + eL_distances_fro(:)), log(1 + lcmv_distances_fro(:)), log(1 + XA_distances_fro(:))};%, log(1 + higgs_distances_fro(:))};
% relative_distances = {log(1 + eL_distances_relative(:)), log(1 + lcmv_distances_relative(:)), log(1 + XA_distances_relative(:))};%, log(1 + higgs_distances_relative(:))};
% corr_distances = {log(1+ eL_distances_corr(:)),log(1+lcmv_distances_corr(:)), log(1+XA_distances_corr(:))};%, log(1+higgs_distances_corr(:))};
% l1_distances = {log(1+eL_distances_l1(:)), log(1+lcmv_distances_l1(:)),log(1+ XA_distances_l1(:))};%, log(1+ higgs_distances_l1(:))};
% 
% % Define colors for the methods (adjust these RGB values as desired)
% colors = [0.2 0.6 0.8;   % Color for eLORETA
%           0.8 0.4 0.2;   % Color for LCMV
%           0.6 0.8 0.2;    % Color for XA
%           0.5 0.5 0.5];  %Color higgs
% 
% % Plot daviolinplot for Frobenius norm distances without scatter data
% % Prepare a new figure for all subplots in a 1x4 layout (excluding Correlation Distance)
% figure;
% 
% % Define font sizes
% titleFontSize = 18;
% labelFontSize = 14;
% legendFontSize = 8;
% 
% % Subplot 1: ROC Curve
% subplot(1, 5, 1);
% plot(eL_X, eL_Y, 'Color', colors(1,:), 'LineWidth', 2); hold on; % eLORETA
% plot(lcmv_X, lcmv_Y, 'Color', colors(2,:), 'LineWidth', 2);      % LCMV
% plot(XA_X, XA_Y, 'Color', colors(3,:), 'LineWidth', 2);          % XA
% %plot(higgs_X, higgs_Y, 'Color', colors(4,:), 'LineWidth', 2);    % Higgs
% legend(['eLORETA (AUC = ' num2str(eL_AUC, '%.2f') ')'], ...
%       ['LCMV (AUC = ' num2str(lcmv_AUC, '%.2f') ')'], ...
%       ['XA (AUC = ' num2str(XA_AUC, '%.2f') ')'], 'Location', 'Best', ...
%       'FontSize', 10, 'FontWeight', 'bold');  % Set legend to bold and size 10
% xlabel('False Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');  % Set axes labels to bold and size 14
% ylabel('True Positive Rate', 'FontSize', 14, 'FontWeight', 'bold');
% title('Combined ROC', 'FontSize', 18, 'FontWeight', 'bold');  % Title to bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks and y-ticks to bold size 14
% hold off;
% 
% % Subplot 2: Frobenius Norm Distance
% subplot(1, 5, 2);
% daviolinplot(fro_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('Frobenius Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('Frobenius Norm', 'FontSize', 18, 'FontWeight', 'bold');       % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 3: Relative Frobenius Error
% subplot(1, 5, 3);
% daviolinplot(relative_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('Frob Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('Frob Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 4: L1 Distance
% subplot(1, 5, 4);
% daviolinplot(l1_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('L1 Distance', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('L1 Norm', 'FontSize', 18, 'FontWeight', 'bold');        % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Subplot 5: Spectral Norm Distance (L2 Distance)
% subplot(1, 5, 5);
% daviolinplot(corr_distances, 'violin', 'half', 'box', 3, 'boxcolors', 'w', ...
%              'xtlabels', methods, 'scatter', 0, ...
%              'violinalpha', 0.7, 'colors', colors, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
%              'boxspacing', 1.2); % Increased boxplot spacing
% ylabel('L1 Relative Error', 'FontSize', 14, 'FontWeight', 'bold');  % Y-label in bold and size 14
% title('L1 Relative Error', 'FontSize', 18, 'FontWeight', 'bold');   % Title in bold and size 18
% grid on;
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set x-ticks (methods) and y-ticks (scale) to bold size 14
% 
% % Adjust the layout for better visualization
% set(gcf, 'Position', [100, 100, 1800, 500]); % Adjust figure size for the 1x5 wide layout
% set(gcf,'Color','w');
% 
% % 
% % %%
% % % Display mean and standard deviation of distances
% % mean_eL_fro = mean(eL_distances_fro);
% % std_eL_fro = std(eL_distances_fro);
% % 
% % mean_lcmv_fro = mean(lcmv_distances_fro);
% % std_lcmv_fro = std(lcmv_distances_fro);
% % 
% % mean_XA_fro = mean(XA_distances_fro);
% % std_XA_fro = std(XA_distances_fro);
% % 
% % fprintf('Frobenius Norm Distances:\n');
% % fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_fro, std_eL_fro);
% % fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_fro, std_lcmv_fro);
% % fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_fro, std_XA_fro);
% % 
% % % Similarly for spectral norm
% % mean_eL_spectral = mean(eL_distances_spectral);
% % std_eL_spectral = std(eL_distances_spectral);
% % 
% % mean_lcmv_spectral = mean(lcmv_distances_spectral);
% % std_lcmv_spectral = std(lcmv_distances_spectral);
% % 
% % mean_XA_spectral = mean(XA_distances_spectral);
% % std_XA_spectral = std(XA_distances_spectral);
% % 
% % fprintf('\nSpectral Norm Distances:\n');
% % fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_spectral, std_eL_spectral);
% % fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_spectral, std_lcmv_spectral);
% % fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_spectral, std_XA_spectral);
% % 
% % % And for relative error
% % mean_eL_relative = mean(eL_distances_relative);
% % std_eL_relative = std(eL_distances_relative);
% % 
% % mean_lcmv_relative = mean(lcmv_distances_relative);
% % std_lcmv_relative = std(lcmv_distances_relative);
% % 
% % mean_XA_relative = mean(XA_distances_relative);
% % std_XA_relative = std(XA_distances_relative);
% % 
% % fprintf('\nRelative Errors:\n');
% % fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_relative, std_eL_relative);
% % fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_relative, std_lcmv_relative);
% % fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_relative, std_XA_relative);
% % 
% % % Display for correlation distance
% % mean_eL_corr = mean(eL_distances_corr);
% % std_eL_corr = std(eL_distances_corr);
% % 
% % mean_lcmv_corr = mean(lcmv_distances_corr);
% % std_lcmv_corr = std(lcmv_distances_corr);
% % 
% % mean_XA_corr = mean(XA_distances_corr);
% % std_XA_corr = std(XA_distances_corr);
% % 
% % fprintf('\nCorrelation Distances:\n');
% % fprintf('eLORETA: Mean = %.4f, Std = %.4f\n', mean_eL_corr, std_eL_corr);
% % fprintf('LCMV: Mean = %.4f, Std = %.4f\n', mean_lcmv_corr, std_lcmv_corr);
% % fprintf('XA: Mean = %.4f, Std = %.4f\n', mean_XA_corr, std_XA_corr);
% 
% 
% %
% % R = voxel_roi_map;
% % preprocessing_velocity;
% % K= parameters.Model.K;
% % L=K;
% % T = parameters.Model.T;
% % K =pinv(K);
% % U_map =K;
% % for j=1:Nw
% %      T_omega(:,:,j) = U_map*T(:,:,j);
% % end
% % parameters.Model.U= T_omega;
% % load('Data/Average_Velocity_ROISpace/GPfit_Delay_Mean.mat');
% % %
% % x_mean = mean(x_avg,2);
% % x_std  = std(x_avg',1)';
% % [e_mean,a_mean,s2_mean]=x2v(x_mean);
% % [e_std,a_std,s2_std]=x2v(x_std);
% % 
% % e_mean =max(e_mean).*ones(size(e_mean));
% % a_mean =max(a_mean).*ones(size(a_mean));
% % s2_mean =max(s2_mean).*ones(size(s2_mean));
% % x_mean = v2x(e_mean,a_mean,s2_mean);
% % 
% % e_std =max(e_std).*ones(size(e_std));
% % a_std =max(a_std).*ones(size(a_std));
% % s2_std =max(s2_std).*ones(size(s2_std));
% % x_std = v2x(e_std,a_mean,s2_std);
% %%
%  S = data_struct.CrossM;
%  parameters.Model = Compact_Model;
% parameters.Compact_Model = Compact_Model;
% parameters.Dimensions = Dimensions;
% Ne = parameters.Dimensions.Ne;
% Nr = parameters.Dimensions.Nr;
% parameters.Dimensions.Nv = Nr;
% Nv = parameters.Dimensions.Nv; % Set voxel to roi
% Nw = parameters.Dimensions.Nw;
% 
% import functions.*
% import functions.auxx.Simulations.*
% import functions.auxx.Simulations.mn_cross.*
% import functions.auxx.DataPreprosessing.*
% S = aveReference(S);
% % S = J;
% %S = Source_PSD;
% freq = data_struct.freqrange(1:47);
% l = zeros(size(S,1),47);
%  Sjj = mn_cross(S, parameters);
% for  i = 1:19
%     for j = 1:47
%         l(i,j) = (diag((Sjj(i,i,j))));
%     end
% end
% plot(freq(1:47),mean(l',2))
% %%
% S = data_struct.CrossM;
% freq = data_struct.freqrange;
% l = zeros(size(S,1), 90);
% 
% % Calculate the l values
% for i = 1:19
%     for j = 1:90
%         l(i,j) = 10 * log10(diag(S(i,i,j)));
%     end
% end
% 
% % Increase smoothness by increasing the number of interpolation points
% x = 1:90; % Original x-values (data points corresponding to the frequency range)
% xq = linspace(1, 90, 500); % Increased number of points for smoother curve
% 
% % Interpolate each row (i.e., each 'i' value) using spline interpolation
% l_smooth = zeros(size(l, 1), 500); % Pre-allocate space for smoothed data
% for i = 1:19
%     l_smooth(i, :) = spline(x, l(i, :), xq); % Interpolate each row
% end
% 
% % Interpolate the frequency range
% freq_smooth = interp1(x, freq(1:90), xq, 'spline'); % Interpolate the frequency range using spline interpolation
% 
% % Plot the smoothed data with increased line width
% figure; % Create a new figure window
% plot(freq_smooth, l_smooth', 'LineWidth', 2) % Set LineWidth to 2 for thicker lines
% xlabel('Frequency (Hz)');
% ylabel('Log Power (dB)');
% title('Smoothed Logarithmic Values of S');
% %%
% S = l;
% [N_c, Nw] = size(l);
% 
% % Define Xi-Omega and Alpha-Omega model functions
% xi_omega_model = @(e1, e2, e3, omega) e1 ./ (1 + e2 * omega.^2).^e3;
% alpha_omega_model = @(a1, a2, a3, a4, omega) a1 ./ (1 + a2 * (omega - a3).^2).^a4;
% 
% % Set up the frequency vector (omega)
% omega = freq; % Adjust the frequency range and number of points (Nw)
% 
% % Initialize storage for fitted parameters
% xi_params = zeros(N_c, 3); % xi_omega parameters for each channel (e1, e2, e3)
% alpha_params = zeros(N_c, 4); % alpha_omega parameters for each channel (a1, a2, a3, a4)
% 
% % Loop over each channel and perform fitting
% for i = 1:N_c
%     % Extract the spectrum for the current channel
%     spectrum = S(i, :); % Spectrum for channel i, of size Nw
% 
%     % Define the combined model function for fitting (7 parameters total)
%     combined_model = @(params, omega) ...
%         xi_omega_model(params(1), params(2), params(3), omega) + ...
%         alpha_omega_model(params(4), params(5), params(6), params(7), omega);
% 
%     % Define the error function (sum of squared errors)
%     % Add penalty for a3 outside the interval [8, 13]
%     error_func = @(params) sum((combined_model(params, omega) - spectrum).^2) + ...
%         1e6 * (params(6) < 8 || params(6) > 13); % Penalize if a3 (params(6)) is outside [8, 13]
% 
%     % Set initial guess for the parameters (7 parameters total)
%     initial_params = [1, 1, 1, 1, 1, 10, 1]; % Example initial guess, set a3 within the range [8, 13]
% 
%     % Transform parameters to ensure positivity (exponential transformation)
%     transform_func = @(x) exp(x); % Exponentiate to ensure positivity
% 
%     % Use fminsearch to minimize the error function
%     options = optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000);
% 
%     % Minimize the error function with transformed parameters
%     fminsearch_error_func = @(transformed_params) error_func(transform_func(transformed_params));
% 
%     fitted_params_transformed = fminsearch(fminsearch_error_func, log(initial_params), options); % Log to apply exp later
% 
%     % Transform the fitted parameters back to positive values
%     fitted_params = transform_func(fitted_params_transformed);
% 
%     % Store the fitted parameters for Xi and Alpha
%     xi_params(i, :) = fitted_params(1:3); % Store Xi-Omega parameters (e1, e2, e3)
%     alpha_params(i, :) = fitted_params(4:7); % Store Alpha-Omega parameters (a1, a2, a3, a4)
% 
%     % Plot the original spectrum and the fitted model
%     figure;
%     plot(omega, spectrum, 'b', 'LineWidth', 1.5); % Original spectrum in blue
%     hold on;
%     plot(omega, combined_model(fitted_params, omega), 'r--', 'LineWidth', 1.5); % Fitted model in red dashed line
%     title(['Channel ' num2str(i) ' Fitting']);
%     xlabel('Frequency (omega)');
%     ylabel('Spectrum Amplitude');
%     legend('Original Spectrum', 'Fitted Model');
%     grid on;
%     hold off;
% end
% 
% % xi_params and alpha_params will now hold the fitted parameters for each channel
% 
% 

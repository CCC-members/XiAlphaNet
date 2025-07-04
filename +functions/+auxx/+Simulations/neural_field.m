


%%
% Import required helper functions
clear all;
close all;
import functions.auxx.ModelVectorization.*;
import guide.Visualization.*;
import functions.auxx.ZeroInflatedModels.*;
import functions.auxx.Refine_Solution.*;
import functions.auxx.OptimizedOperations.*;
Cortex = load("templates/Cortex.mat");

% Path to the JSON file with model result metadata
json_path = '/mnt/Store/Ronaldo/dev/Data/NewFolder/XIALPHANET.json';
[dataset_dir, ~, ~] = fileparts(json_path);
dataset = jsondecode(fileread(json_path));
dataset.Location = dataset_dir;
parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));

parameters.Model = parameters.Compact_Model;
parameters.Compact_Model = parameters.Compact_Model;
parameters.Dimensions = parameters.Dimensions;
parameters.Dimensions.Nv = parameters.Dimensions.Nr;
% Clear temporary variables to save memory
clear Compact_Model Model Dimensions;

% Import additional functions
import functions.*;
import functions.auxx.*;
import functions.auxx.BayesOptimization.*;
import functions.auxx.CrossSpectrum.*;
import functions.auxx.ModelVectorization.*;
import functions.auxx.Simulations.*;
import functions.auxx.TOperator.*;
import tools.*;
import functions.auxx.DataPreprocessing.*;
import functions.auxx.Simulations.inverse.*;
import functions.auxx.Simulations.private.*;

% Set up simulation parameters
Ne = parameters.Dimensions.Ne;  % Number of electrodes
Nr = parameters.Dimensions.Nr;  % Number of ROIs
Nv = parameters.Dimensions.Nr;  % Voxel dimension = ROI dimension
Nw = parameters.Dimensions.Nw;  % Number of frequency bins
Nsub = 100;                      % Number of simulations
N_wishart = 1;
conn_spec = norm(parameters.Model.C, 'fro');

disp("--> Estimating source cross-spectrum");

% Xi-AlphaNET properties
properties.model_params.nFreqs = Nw;
properties.model_params.BayesIter_Delay = 50;
properties.model_params.BayesIter_Reg1 = 20;
properties.model_params.BayesIter_Reg2 = 100;
properties.model_params.Nrand1 = 10;
properties.model_params.Nrand2 = 50;
properties.model_params.delay.lambda_space_cd = [[0.4, 1.6]; [10^(-10), 1/conn_spec]];
properties.general_params.parallel.conn_delay = 1;

properties.model_params.stoch1 = 1;
properties.model_params.stoch2 = 1;
properties.model_params.tensor_field.default = 0;

% Load model and transformation matrices
L = parameters.Model.K;  % Transformation matrix for cross-spectrum
K = parameters.Model.K;
R = parameters.Compact_Model.R;

% Select subjects for simulation
dir_data = '/mnt/Store/Ronaldo/dev/Data/norms';
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
selected_folders = subject_folders(randperm(length(subject_folders), Nsub));

% Define frequency bands
bands = {
    'Delta', 1, 4;
    'Theta', 4, 8;
    'Alpha', 8, 13;
    'Beta', 13, 30;
    'Gamma', 30, 50
    };

Nbands = size(bands, 1);

% Pre-allocate variables
Nw = 47; % Assuming freq is defined elsewhere
Nroi = 360;        % Adjust if needed

xl = zeros(Nroi, Nw, Nsub);
l  = zeros(Nroi, Nw, Nsub);
mf = zeros(Nroi, Nw, Nsub);
% Frobenius and L1 norm benchmark errors for each subject
XA_Spectrum_Fro = zeros(Nsub, 1);
XA_Spectrum_FroR = zeros(Nsub, 1);
XA_Spectrum_L1 = zeros(Nsub, 1);
XA_Spectrum_L1R = zeros(Nsub, 1);

mn_Spectrum_Fro = zeros(Nsub, 1);
mn_Spectrum_FroR = zeros(Nsub, 1);
mn_Spectrum_L1 = zeros(Nsub, 1);
mn_Spectrum_L1R = zeros(Nsub, 1);

XA_Coherence_Fro = zeros(Nsub, 1);
XA_Coherence_FroR = zeros(Nsub, 1);
XA_Coherence_L1 = zeros(Nsub, 1);
XA_Coherence_L1R = zeros(Nsub, 1);

mn_Coherence_Fro = zeros(Nsub, 1);
mn_Coherence_FroR = zeros(Nsub, 1);
mn_Coherence_L1 = zeros(Nsub, 1);
mn_Coherence_L1R = zeros(Nsub, 1);

XA_Phase_Fro = zeros(Nsub, 1);
XA_Phase_FroR = zeros(Nsub, 1);
XA_Phase_L1 = zeros(Nsub, 1);
XA_Phase_L1R = zeros(Nsub, 1);

mn_Phase_Fro = zeros(Nsub, 1);
mn_Phase_FroR = zeros(Nsub, 1);
mn_Phase_L1 = zeros(Nsub, 1);
mn_Phase_L1R = zeros(Nsub, 1);


% Simulation loop
for j = 1:Nsub
    tic;
    % Load subject data
    subject_folder = selected_folders(j).name;
    mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
    data_struct = load(mat_file_path);

    freq = data_struct.data_struct.freqrange(1:Nw);
    parameters.Data.freq = freq;

    % Source cross-spectrum simulation
    %%
    disp('->> Estimating Source Cross Using a Neural Mass Simulation');
    [Sjj, ~] = functions.auxx.Simulations.neural_mass_simulation(1);

    % Generate scalp-level cross-spectrum with noise

    N_wishart = 4;
    disp('->> Simulating Scalp Cross Wishart Noise + FModel');
    parfor i = 1:Nw
        Sjj_cross(:, :, i) = generate_complex_wishart(Sjj(:, :, i), N_wishart);
        Svv_cross(:, :, i) = L * Sjj_cross(:, :, i) * L';  % Apply transformation
    end
    [Log_Spec,freq]= log_spectrum(Svv_cross,freq);
    plot(freq,Log_Spec')

    %%
    toc;

    % Xi-AlphaNET estimation
    data.Cross = Svv_cross;
    data.age = 25;  % Subject age
    data.freq = freq;

    disp('->> Xi-AlphaNeT Inverse Solution');
    [x, ~, G, x0] = Xi_ALphaNET(properties, data, parameters);
    source_act_cross = functions.auxx.CrossSpectrum.eval_source_conn(x.Solution, data.freq, parameters.Model.R, properties, parameters);
    XA_Sjj_cross = source_act_cross.Cross.Full;

    % Additional processing (mean cross-spectrum)
    % Calculate the mean cross-spectrum across frequencies
    mn_Sjj_cross = mn_cross(Svv_cross, K, 0);

    disp('->> eLORETA Processing...');

    % Compute eLORETA source cross-spectra for each frequency bin
    parfor i = 1:Nw
        source = inverse(Svv_cross(:,:,i), L);
        eL_Sjj_cross(:,:,i) = source.eloreata.Sjj;
    end

    % Initialize mean power array for each ROI and frequency
    mn_power = zeros(Nroi, Nw);

    % Extract and log-transform diagonal power for each ROI
    for i = 1:Nw
        xl(:, i, j)= log10(real(diag(XA_Sjj_cross(:, :, i))));
        % Original cross-spectrum power (log-transformed)
        l(:, i, j) = log10(real(diag(Sjj_cross(:, :, i))));
        % eLORETA power (not log-transformed)
        lo(:, i) = real(diag(eL_Sjj_cross(:, :, i)));
        % Mean cross-spectrum power
        mn_power(:, i) = real(diag(mn_Sjj_cross(:, :, i)));
    end

    % Fit FOOOF and xi-alpha models per ROI
    parfor roi = 1:Nroi
        roi
        % Fit the FOOOF model to the mean power spectrum
        mn_fooof = fooof_matlab(freq, mn_power(roi, :), [min(freq), max(freq)], 1);
        mf(roi, :, j) = mn_fooof.model_fit(:).';

        % Fit the xi-alpha model to eLORETA diagonal power
        [params_xa_lo, fit_xa_lo] = functions.auxx.GenerateSourceSample.fit_xi_alpha_multi(lo(roi,:), freq,10, 0);
        el_fit_xa(roi,:,j) = fit_xa_lo;
    end

    % Store ground-truth cross-spectrum (used for comparisons)
    ground_truth_Sjj = Sjj_cross;

    % ----- Spectral Comparison -----
    % Compare estimated cross-spectra to ground-truth using Frobenius and L1 norms

    % Notation:
    % - "XA" corresponds to your Xi-AlphaNET or a different model; rename as needed
    % - "mn" corresponds to the mean cross-spectrum
    % - "eL" corresponds to the eLORETA estimates

    % Frobenius Norm Comparison
    XA_Spectrum_Fro(j)   = tensor_norm(XA_Sjj_cross - ground_truth_Sjj, 2);
    mn_Spectrum_Fro(j)   = tensor_norm(mn_Sjj_cross - ground_truth_Sjj, 2);
    eL_Spectrum_Fro(j)   = tensor_norm(eL_Sjj_cross - ground_truth_Sjj, 2);

    % Relative Frobenius Norm (normalized by ground-truth)
    XA_Spectrum_FroR(j)  = XA_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    mn_Spectrum_FroR(j)  = mn_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    eL_Spectrum_FroR(j)  = eL_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);

    % L1 Norm Comparison
    XA_Spectrum_L1(j)    = tensor_norm(XA_Sjj_cross - ground_truth_Sjj, 1);
    mn_Spectrum_L1(j)    = tensor_norm(mn_Sjj_cross - ground_truth_Sjj, 1);
    eL_Spectrum_L1(j)    = tensor_norm(eL_Sjj_cross - ground_truth_Sjj, 1);

    % Relative L1 Norm (normalized by ground-truth)
    XA_Spectrum_L1R(j)   = XA_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    mn_Spectrum_L1R(j)   = mn_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    eL_Spectrum_L1R(j)   = eL_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);

    % ----- Coherence Comparison -----
    % Compare coherence matrices for each model
    coherence_Sjj = coherence(ground_truth_Sjj);
    coherence_XA  = coherence(XA_Sjj_cross);
    coherence_mn  = coherence(mn_Sjj_cross);
    coherence_eL  = coherence(eL_Sjj_cross);

    % Frobenius Norm Comparison for Coherence
    XA_Coherence_Fro(j)  = tensor_norm(coherence_XA - coherence_Sjj, 2);
    mn_Coherence_Fro(j)  = tensor_norm(coherence_mn - coherence_Sjj, 2);
    eL_Coherence_Fro(j)  = tensor_norm(coherence_eL - coherence_Sjj, 2);

    % Relative Frobenius Norm
    XA_Coherence_FroR(j) = XA_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    mn_Coherence_FroR(j) = mn_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    eL_Coherence_FroR(j) = eL_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);

    % L1 Norm Comparison for Coherence
    XA_Coherence_L1(j)   = tensor_norm(coherence_XA - coherence_Sjj, 1);
    mn_Coherence_L1(j)   = tensor_norm(coherence_mn - coherence_Sjj, 1);
    eL_Coherence_L1(j)   = tensor_norm(coherence_eL - coherence_Sjj, 1);

    % Relative L1 Norm
    XA_Coherence_L1R(j)  = XA_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    mn_Coherence_L1R(j)  = mn_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    eL_Coherence_L1R(j)  = eL_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);

    % ----- Phase Comparison -----
    % Compare phase angles across models
    angle_Sjj = angle(ground_truth_Sjj);
    angle_XA  = angle(XA_Sjj_cross);
    angle_mn  = angle(mn_Sjj_cross);
    angle_eL  = angle(eL_Sjj_cross);

    % Frobenius Norm Comparison for Phase
    XA_Phase_Fro(j)      = tensor_norm(angle_XA - angle_Sjj, 2);
    mn_Phase_Fro(j)      = tensor_norm(angle_mn - angle_Sjj, 2);
    eL_Phase_Fro(j)      = tensor_norm(angle_eL - angle_Sjj, 2);

    % Relative Frobenius Norm
    XA_Phase_FroR(j)     = XA_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    mn_Phase_FroR(j)     = mn_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    eL_Phase_FroR(j)     = eL_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);

    % L1 Norm Comparison for Phase
    XA_Phase_L1(j)       = tensor_norm(angle_XA - angle_Sjj, 1);
    mn_Phase_L1(j)       = tensor_norm(angle_mn - angle_Sjj, 1);
    eL_Phase_L1(j)       = tensor_norm(angle_eL - angle_Sjj, 1);

    % Relative L1 Norm
    XA_Phase_L1R(j)      = XA_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    mn_Phase_L1R(j)      = mn_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    eL_Phase_L1R(j)      = eL_Phase_L1(j) / tensor_norm(angle_Sjj, 1);


end

% -----------------------------------
%% Colors
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*
colors = [0.2 0.6 0.8;   % eLORETA
    0.8 0.4 0.2;   % LCMV
    0.6 0.8 0.2;   % Xi-AlphaNET
    0.5 0.5 0.5];  % gray

color_xa = colors(3, :);
color_foof = colors(1, :);
color_elxa = colors(2,:);
% Overall MSE
MSE_xl = zeros(Nsub, Nroi);
MSE_mf = zeros(Nsub, Nroi);
MSE_lo = zeros(Nsub, Nroi);

for j = 1:Nsub
    for roi = 1:Nroi
        l_roi  = squeeze(l(roi,:,j));
        xl_roi = squeeze(xl(roi,:,j));
        mf_roi = squeeze(mf(roi,:,j));
        lo_roi = squeeze(el_fit_xa(roi,:,j));

        MSE_xl(j, roi) = mean((xl_roi - l_roi).^2, 'omitnan');
        MSE_mf(j, roi) = mean((mf_roi - l_roi).^2, 'omitnan');
        MSE_lo(j, roi) = mean((lo_roi - l_roi).^2, 'omitnan');
    end
end

mean_MSE_xl = mean(MSE_xl, 2, 'omitnan');
mean_MSE_mf = mean(MSE_mf, 2, 'omitnan');
mean_MSE_lo = mean(MSE_lo, 2, 'omitnan');

% Paired t-test
[~, p_ttest] = ttest(mean_MSE_xl, mean_MSE_mf);
disp(['Paired t-test p-value (MSE comparison): ', num2str(p_ttest)]);

% Band-specific MSEs
freq_band_indices = cell(Nbands,1);
for b = 1:Nbands
    freq_band_indices{b} = find(freq >= bands{b,2} & freq <= bands{b,3});
end

MSE_xl_band = zeros(Nsub, Nroi, Nbands);
MSE_mf_band = zeros(Nsub, Nroi, Nbands);
MSE_eL_band = zeros(Nsub, Nroi, Nbands);

for j = 1:Nsub
    for roi = 1:Nroi
        l_roi  = squeeze(l(roi,:,j));
        xl_roi = squeeze(xl(roi,:,j));
        mf_roi = squeeze(mf(roi,:,j));
        lo_roi = squeeze(el_fit_xa(roi,:,j));
        for b = 1:Nbands
            idx = freq_band_indices{b};
            MSE_xl_band(j, roi, b) = mean((xl_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_mf_band(j, roi, b) = mean((mf_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_eL_band(j, roi, b) = mean((lo_roi(idx)-l_roi(idx)).^2, 'omitnan');
        end
    end
end

mean_MSE_xl_band = squeeze(mean(MSE_xl_band,2,'omitnan'));
mean_MSE_mf_band = squeeze(mean(MSE_mf_band,2,'omitnan'));
mean_MSE_eL_band = squeeze(mean(MSE_eL_band,2,'omitnan'));


% Overall violin plot
errors_combined = {mean_MSE_xl, mean_MSE_mf,mean_MSE_lo};
methods = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha'};
colors_overall = [color_xa; color_foof;color_elxa];

figure('Color','w');
daviolinplot(errors_combined, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_overall, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('Mean Squared Error','FontSize',14,'FontWeight','bold');
title('Overall Reconstruction Errors','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);

% Plot band-specific violin plots with stars for minimum MSE

selected_bands = 1:(Nbands-1); % Exclude Gamma
errors_band = {};
xtlabels_band = {};
colors_band = [];

for idx = 1:numel(selected_bands)
    b = selected_bands(idx);

    % Prepare cell array for daviolinplot
    errors_band{end+1} = log10(mean_MSE_xl_band(:,b));
    errors_band{end+1} = log10(mean_MSE_mf_band(:,b));
    errors_band{end+1} = log10(mean_MSE_eL_band(:,b));        
    band_name = bands{b,1};
    if iscell(band_name), band_name = band_name{1}; end
    if isstring(band_name), band_name = char(band_name); end

    xtlabels_band{end+1} = [band_name, ' (\xi-\alphaNET)'];
    xtlabels_band{end+1} = [band_name, ' (MNE+FOOOF)'];
    xtlabels_band{end+1} = [band_name, ' (eLORETA+\xi-\alpha)'];

    colors_band = [colors_band; color_xa; color_foof;color_elxa];
end

figure('Color','w');
daviolinplot(errors_band, 'violin', 'half', 'box', 3, ...
    'xtlabels', xtlabels_band, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_band, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('Mean Squared Error','FontSize',14,'FontWeight','bold');
title('Band-Specific Reconstruction Errors (Gamma Excluded)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);
xtickangle(45);
%%
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*

% Define colors
colors = [0.2 0.6 0.8;   % eLORETA
          0.8 0.4 0.2;   % LCMV
          0.6 0.8 0.2;]; % Xi-AlphaNET

color_xa = colors(3, :);
color_foof = colors(1, :);
color_elxa = colors(2,:);

% ----- Coherence Performance -----
% Average across ROI if necessary and apply log transform
if size(XA_Coherence_FroR,2) > 1
    XA_Coherence_FroR_mean = mean(XA_Coherence_FroR, 2, 'omitnan');
else
    XA_Coherence_FroR_mean = XA_Coherence_FroR;
end
XA_Coherence_FroR_log = log10(XA_Coherence_FroR_mean);

if size(mn_Coherence_FroR,2) > 1
    mn_Coherence_FroR_mean = mean(mn_Coherence_FroR, 2, 'omitnan');
else
    mn_Coherence_FroR_mean = mn_Coherence_FroR;
end
mn_Coherence_FroR_log = log10(mn_Coherence_FroR_mean);

if size(eL_Coherence_FroR,2) > 1
    eL_Coherence_FroR_mean = mean(eL_Coherence_FroR, 2, 'omitnan');
else
    eL_Coherence_FroR_mean = eL_Coherence_FroR;
end
eL_Coherence_FroR_log = log10(eL_Coherence_FroR_mean);

errors_coherence = {XA_Coherence_FroR_log, mn_Coherence_FroR_log, eL_Coherence_FroR_log};
methods_coherence = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha'};
colors_coherence = [color_xa; color_foof; color_elxa];

figure('Color','w');
daviolinplot(errors_coherence, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods_coherence, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_coherence, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('log_{10}(Relative Frobenius Norm)','FontSize',14,'FontWeight','bold');
title('Performance Comparison (Coherence)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);

% ----- Phase Performance -----
if size(XA_Phase_FroR,2) > 1
    XA_Phase_FroR_mean = mean(XA_Phase_FroR, 2, 'omitnan');
else
    XA_Phase_FroR_mean = XA_Phase_FroR;
end
XA_Phase_FroR_log = log10(XA_Phase_FroR_mean);

if size(mn_Phase_FroR,2) > 1
    mn_Phase_FroR_mean = mean(mn_Phase_FroR, 2, 'omitnan');
else
    mn_Phase_FroR_mean = mn_Phase_FroR;
end
mn_Phase_FroR_log = log10(mn_Phase_FroR_mean);

if size(eL_Phase_FroR,2) > 1
    eL_Phase_FroR_mean = mean(eL_Phase_FroR, 2, 'omitnan');
else
    eL_Phase_FroR_mean = eL_Phase_FroR;
end
eL_Phase_FroR_log = log10(eL_Phase_FroR_mean);

errors_phase = {XA_Phase_FroR_log, mn_Phase_FroR_log, eL_Phase_FroR_log};
methods_phase = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha'};
colors_phase = [color_xa; color_foof; color_elxa];

figure('Color','w');
daviolinplot(errors_phase, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods_phase, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_phase, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('log_{10}(Relative Frobenius Norm)','FontSize',14,'FontWeight','bold');
title('Performance Comparison (Phase)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);

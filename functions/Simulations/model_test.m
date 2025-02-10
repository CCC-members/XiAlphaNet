%Povilas Karvelis (2024). daviolinplot - beautiful violin and raincloud plots 
% (https://github.com/povilaskarvelis/DataViz/releases/tag/v3.2.4), GitHub. Retrieved September 20, 2024.


clear all 
clc;
load("Data\Model_Parameters\parameters.mat")
load('Data\Model_Parameters\x_avg.mat');
% Setting initialiation of te simulations  
Ne=19;
Nr= 360;
Nf = 49;
Nsim = 1;
%
Svv_cross=zeros(Ne,Ne,Nf,Nsim);  % structure to safe the cross 
Sjj_cross=zeros(Nr,Nr,Nf,Nsim);  % structure to safe the cross 
eL_Sjj_cross= zeros(Nr,Nr,Nf,Nsim);
lcmv_Sjj_cross= zeros(Nr,Nr,Nf,Nsim);
XA_Sjj_cross= zeros(Nr,Nr,Nf,Nsim);
higgs_Sjj_cross = zeros(Nr,Nr,Nf,Nsim);
%
R = voxel_roi_map;
preprocessing_velocity;
K= parameters.Model.K;
L=K;
T = parameters.Model.T;
K =pinv(K);
U_map =K;
for j=1:Nf
     T_omega(:,:,j) = U_map*T(:,:,j);
end
parameters.Model.U= T_omega;
load('Data/Average_Velocity_ROISpace/GPfit_Delay_Mean.mat');
%
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
%

% Generate a random sample of plausible crospectrum in the scalp 
Svv = parameters.Data.Cross;
Sjj = mn_cross(Svv,parameters);
N_wishart = 100;
for j=1:Nsim 
    for i=1:Nf
    Sjj_cross(:,:,i,j) = generate_complex_wishart(Sjj(:,:,i), N_wishart);
    Svv_cross(:,:,i,j) = L*Sjj_cross(:,:,i,j)*L';
    end
end

% 

for j=1:Nsim 
    j
    tic
    parameters.Data.Cross = Svv_cross(:,:,:,j);
    [S,~] =  Xi_ALphaNET(parameters);
    XA_Sjj_cross(:,:,:,j) = S;
    for i=1:8
        i
        % Find the etimation using eLoreta and LCMV 
        source = inverse(Svv_cross(:,:,i,j),L);
        eL_Sjj_cross(:,:,i,j) = source.eloreata.Sjj;
        lcmv_Sjj_cross(:,:,i,j) = source.lcmv.Sjj;
        higgs_Sjj_cross(:,:,i,j) = call_higgs(Svv_cross(:,:,i,j),L);
    end
    toc
end
%%



%% Read Data 
 Svv_cross = data_struct.Svv_cross;
Sjj_cross = data_struct.Sjj_cross;
eL_Sjj_cross = data_struct.eL_Sjj_cross;
lcmv_Sjj_cross = data_struct.lcmv_Sjj_cross;
XA_Sjj_cross = data_struct.XA_Sjj_cross;
data_struct2.Svv_cross=Svv_cross;
data_struct2.Sjj_cross=Sjj_cross;
data_struct2.eL_Sjj_cross=eL_Sjj_cross;
data_struct2.lcmv_Sjj_cross=lcmv_Sjj_cross;
data_struct2.XA_Sjj_cross=XA_Sjj_cross;

% %%
% Define your parameters
Nr = size(Sjj_cross, 1); % Number of regions
Nf = 25;%size(Sjj_cross, 3); % Number of frequencies
%Nsim = 50; % Number of simulations

% Precompute the upper triangular indices (excluding the diagonal)
idx_triu = triu(true(Nr), 1); % Only needs to be calculated once

% Preallocate memory for scores and distances based on estimated sizes
num_elements = sum(idx_triu(:)); % Number of upper triangular elements
total_size = Nsim * Nf * num_elements;

all_labels = zeros(total_size, 1);
all_eL_scores = zeros(total_size, 1);
all_lcmv_scores = zeros(total_size, 1);
all_XA_scores = zeros(total_size, 1);
all_higgs_scores = zeros(total_size, 1);


% Preallocate distances
eL_distances_fro = zeros(Nsim * Nf, 1);
lcmv_distances_fro = zeros(Nsim * Nf, 1);
XA_distances_fro = zeros(Nsim * Nf, 1);
higgs_distances_fro = zeros(Nsim * Nf, 1);


eL_distances_relative = zeros(Nsim * Nf, 1);
lcmv_distances_relative = zeros(Nsim * Nf, 1);
XA_distances_relative = zeros(Nsim * Nf, 1);
higgs_distances_relative = zeros(Nsim * Nf, 1);


eL_distances_corr = zeros(Nsim * Nf, 1);
lcmv_distances_corr = zeros(Nsim * Nf, 1);
XA_distances_corr = zeros(Nsim * Nf, 1);
higgs_distances_corr = zeros(Nsim * Nf, 1);

eL_distances_l1 = zeros(Nsim * Nf, 1);
lcmv_distances_l1 = zeros(Nsim * Nf, 1);
XA_distances_l1 = zeros(Nsim * Nf, 1);
higgs_distances_l1 = zeros(Nsim * Nf, 1);

% Counters for storing data efficiently
score_counter = 1;
distance_counter = 1;
%%
for n = 1:1
    disp(n); % Display the simulation number
    for f = 1:Nf
        % Extract true and estimated covariance matrices (pre-compute abs)
        Sjj_true = abs((Sjj_cross(:, :, f, n)));
        eL_Sjj_est = abs((eL_Sjj_cross(:, :, f, n)));
        lcmv_Sjj_est = abs((lcmv_Sjj_cross(:, :, f, n)));
        XA_Sjj_est = abs((XA_Sjj_cross(:, :, f, n)));
        higgs_Sjj_est = abs((higgs_Sjj_cross(:, :, f, n)));

        % Flatten the upper triangular parts of the matrices
        Sjj_true_vec = (Sjj_true(idx_triu));
        eL_Sjj_est_vec = (eL_Sjj_est(idx_triu));
        lcmv_Sjj_est_vec = (lcmv_Sjj_est(idx_triu));
        XA_Sjj_est_vec = (XA_Sjj_est(idx_triu));
        higgs_Sjj_est_vec = (higgs_Sjj_est(idx_triu));

        % Threshold for ground truth labels
        threshold = mean(Sjj_true_vec);
        labels = Sjj_true_vec > threshold;

        % Store labels and scores
        num_elements = numel(Sjj_true_vec);
        all_labels(score_counter:score_counter + num_elements - 1) = labels;
        all_eL_scores(score_counter:score_counter + num_elements - 1) = eL_Sjj_est_vec;
        all_lcmv_scores(score_counter:score_counter + num_elements - 1) = lcmv_Sjj_est_vec;
        all_XA_scores(score_counter:score_counter + num_elements - 1) = XA_Sjj_est_vec;
        all_higgs_scores(score_counter:score_counter + num_elements - 1) = higgs_Sjj_est_vec;
        score_counter = score_counter + num_elements;

        % Compute Frobenius norm and relative error
        eL_distances_fro(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro');
        lcmv_distances_fro(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro');
        XA_distances_fro(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro');
        higgs_distances_fro(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro');

        eL_distances_relative(distance_counter) = norm(Sjj_true - eL_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        lcmv_distances_relative(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        XA_distances_relative(distance_counter) = norm(Sjj_true - XA_Sjj_est, 'fro') / norm(Sjj_true, 'fro');
        higgs_distances_relative(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 'fro') / norm(Sjj_true, 'fro');

        % Compute correlation distances
        eL_distances_corr(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1) / norm(Sjj_true, 1);
        lcmv_distances_corr(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1) / norm(Sjj_true, 1);
        XA_distances_corr(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1) / norm(Sjj_true, 1);
        higgs_distances_corr(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1) / norm(Sjj_true, 1);
        % L1 norm
        eL_distances_l1(distance_counter) = norm(Sjj_true - eL_Sjj_est, 1);
        lcmv_distances_l1(distance_counter) = norm(Sjj_true - lcmv_Sjj_est, 1);
        XA_distances_1(distance_counter) = norm(Sjj_true - XA_Sjj_est, 1);
        higgs_distances_1(distance_counter) = norm(Sjj_true - higgs_Sjj_est, 1);
        distance_counter = distance_counter + 1;
    end
end

%%
% Compute ROC curves using all data
[eL_X, eL_Y, ~, eL_AUC] = perfcurve(all_labels, all_eL_scores, true);
[lcmv_X, lcmv_Y, ~, lcmv_AUC] = perfcurve(all_labels, all_lcmv_scores, true);
[XA_X, XA_Y, ~, XA_AUC] = perfcurve(all_labels, all_XA_scores, true);
[higgs_X, higgs_Y, ~, higgs_AUC] = perfcurve(all_labels, all_higgs_scores, true);



%%
%here

% Prepare data for daviolinplot of distances
methods = {'eLORETA', 'LCMV', 'XA', 'higgs'};

% Ensure data are column vectors and apply the log(1 + x) transformation
fro_distances = {log(1 + eL_distances_fro(:)), log(1 + lcmv_distances_fro(:)), log(1 + XA_distances_fro(:)), log(1 + higgs_distances_fro(:))};
relative_distances = {log(1 + eL_distances_relative(:)), log(1 + lcmv_distances_relative(:)), log(1 + XA_distances_relative(:)), log(1 + higgs_distances_relative(:))};
corr_distances = {log(1+ eL_distances_corr(:)),log(1+lcmv_distances_corr(:)), log(1+XA_distances_corr(:)), log(1+higgs_distances_corr(:))};
l1_distances = {log(1+eL_distances_l1(:)), log(1+lcmv_distances_l1(:)),log(1+ XA_distances_l1(:)), log(1+ higgs_distances_l1(:))};

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
plot(higgs_X, higgs_Y, 'Color', colors(4,:), 'LineWidth', 2);    % Higgs
legend(['eLORETA (AUC = ' num2str(eL_AUC, '%.2f') ')'], ...
      ['LCMV (AUC = ' num2str(lcmv_AUC, '%.2f') ')'], ...
      ['XA (AUC = ' num2str(XA_AUC, '%.2f') ')'],['higgs (AUC = ' num2str(higgs_AUC, '%.2f') ')'] , 'Location', 'Best', ...
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

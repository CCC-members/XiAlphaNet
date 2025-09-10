clear; clc;

% Directory containing .mat files results
json_path = '/mnt/Develop/Ronaldo/dev/test262025_delete/XIALPHANET.json';
[dataset_dir, ~, ~] = fileparts(json_path);
dataset = jsondecode(fileread(json_path));
dataset.Location = dataset_dir;

import templates.*
load("templates/mylin_data.mat") % myelin_data.Age, .myelin, .upper, .lower

% === Extract conduction delays and ages ===
delays = []; ages = [];
index = 1;
for i = 1:length(dataset.Participants)
    participant = dataset.Participants(i);
    if isequal(participant.Status, 'Completed')
        Part_Info = jsondecode(fileread(fullfile(dataset.Location, participant.SubID, participant.FileInfo)));
        D = load(fullfile(dataset.Location, participant.SubID, Part_Info.Delay_Matrix));
        delays(index) = 1000 * mean(D.Delay_Matrix(:));
        ages(index) = participant.Age;
        index = index + 1;
    end
end

% === Clean and preprocess ===
delays = delays(:); ages = ages(:);
valid = ~isnan(delays) & ~isnan(ages);
delays = delays(valid); ages = ages(valid);
%plot(ages,delays,"*")
% Z-score outlier removal
z = abs(zscore(delays));
delays = delays(z < 3);
ages = ages(z < 3);

[ages, sortIdx] = sort(ages);
delays = delays(sortIdx);
% === FIGURE 1: Conduction Delay ===
log_delays = delays;
X = [ones(size(ages)), ages, ages.^2];
[b1, stats1] = robustfit(X(:,2:end), log_delays);
ages_fit = linspace(min(ages), max(ages), 200)';
X_fit = [ones(size(ages_fit)), ages_fit, ages_fit.^2];
log_y_fit = X_fit * b1;
y_fit = log_y_fit;
se_log = sqrt(sum((X_fit * stats1.covb) .* X_fit, 2));
upper = log_y_fit + se_log;
lower = log_y_fit - se_log;

figure('Color','w','Units','normalized','Position',[0.2 0.2 0.6 0.6]);

ax_main = axes('Position', [0.15 0.15 0.65 0.62]); hold on;
h_fill = fill([ages_fit; flipud(ages_fit)], [upper; flipud(lower)], 'r', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
h_fit = plot(ages_fit, y_fit, 'r--', 'LineWidth', 2);
xlabel('Age (Years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conduction Delay (ms)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_main, 'FontSize', 13, 'FontWeight', 'bold'); xlim([5 95]); grid on;

% KDE (top)
ax_top = axes('Position', [0.15 0.79 0.65 0.12]);
[fx, xgrid] = ksdensity(ages);
fill(ax_top, xgrid, fx, [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'k'); hold on;
mu_age = mean(ages); sigma_age = std(ages);
gauss_top = normpdf(xgrid, mu_age, sigma_age);
plot(ax_top, xgrid, gauss_top, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
axis(ax_top, 'tight'); xlim([5 95]);
set(ax_top, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

% KDE (right)
ax_right = axes('Position', [0.82 0.15 0.12 0.62]);
[fy, ygrid] = ksdensity(delays);
fill(ax_right, fy, ygrid, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'r'); hold on;
mu_delay = mean(delays); sigma_delay = std(delays);
h_gauss = plot(ax_right, normpdf(ygrid, mu_delay, sigma_delay), ygrid, '--', ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
axis(ax_right, 'tight');
set(ax_right, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

legend(ax_main, [h_fit, h_fill, h_gauss], ...
    {'Delay Fit', 'Estimator Uncertainty', 'Gaussian Fit (Delay)'}, ...
    'Location', 'southwest', 'FontSize', 12);

annotation('textbox', [0.15, 0.93, 0.7, 0.05], ...
    'String', sprintf('Conduction Delay vs Age (p = %.3g)', stats1.p(3)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');

% === FIGURE 2: Myelin Estimate ===
log_myelin = log(1./delays.^2);
pos = ages < 83;
ages = ages(pos);
log_myelin = log_myelin(pos);
X = [ones(size(ages)), ages, ages.^2];
ages_fit = linspace(0.3, 83, 100)';
X_fit = [ones(size(ages_fit)), ages_fit, ages_fit.^2];
[b2, stats2] = robustfit(X(:,2:end), log_myelin);
log_y_fit = X_fit * b2;
y_fit = exp(log_y_fit);
se_log = sqrt(sum((X_fit * stats2.covb) .* X_fit, 2));
upper = exp(log_y_fit + se_log);
lower = exp(log_y_fit - se_log);

Z = max(upper) - min(lower);
y_norm = (y_fit - min(lower)) / Z;
upper = (upper - min(lower)) / Z;
lower = (lower - min(lower)) / Z;

figure('Color','w','Units','normalized','Position',[0.2 0.2 0.6 0.6]);

ax_main = axes('Position', [0.15 0.15 0.65 0.62]); hold on;
h_fill = fill([ages_fit; flipud(ages_fit)], [upper; flipud(lower)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h_fit = plot(ages_fit, y_norm, 'r--', 'LineWidth', 2);
h_data = plot(myelin_data.Age, myelin_data.myelin, 'b--', 'LineWidth', 2);
h_unc = fill([myelin_data.Age, fliplr(myelin_data.Age)], ...
     [myelin_data.upper, fliplr(myelin_data.lower)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Age (Years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Estimated Myelin (1/\tau^2)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_main, 'FontSize', 13, 'FontWeight', 'bold'); xlim([0.3 83]); grid on;

% KDE (top)
ax_top = axes('Position', [0.15 0.79 0.65 0.12]);
[fx, xgrid] = ksdensity(ages);
fill(ax_top, xgrid, fx, [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'k'); hold on;
mu_age = mean(ages); sigma_age = std(ages);
gauss_top = normpdf(xgrid, mu_age, sigma_age);
plot(ax_top, xgrid, gauss_top, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
axis(ax_top, 'tight'); xlim([0.3 83]);
set(ax_top, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

% KDE (right)
ax_right = axes('Position', [0.82 0.15 0.12 0.62]);
[fy, ygrid] = ksdensity(y_norm);
fill(ax_right, fy, ygrid, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'r'); hold on;
mu_myelin = mean(y_norm); sigma_myelin = std(y_norm);
h_gauss = plot(ax_right, normpdf(ygrid, mu_myelin, sigma_myelin), ygrid, '--', ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
axis(ax_right, 'tight'); ylim([0 1]);
set(ax_right, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

legend(ax_main, [h_fit, h_fill, h_data, h_unc, h_gauss], ...
    {'Estimated Myelin', 'Estimator Uncertainty', ...
     'Myelin Data', 'Myelin Uncertainty', 'Gaussian Fit'}, ...
    'Location', 'southwest', 'FontSize', 12);

annotation('textbox', [0.15, 0.93, 0.7, 0.05], ...
    'String', sprintf('Estimated Myelin vs Myelin Data (p = %.5g)', stats2.p(3)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');


%% Spectral components vs delays
clc; clear all;

% Import required helper functions
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*

% Define parameters for analysis
prc = 90;             % Percentile threshold (e.g., for visualizing distributions)
cross_index = 0;      % Cross-validation index (0 = not using CV)
num_groups = 1;       % Number of age groups (for age stratification)
mode = 0;             % 0 = amplitude plots, 1 = zero-inflation probability plots
age_min = 0;          % Minimum age for inclusion
age_max = 100;        % Maximum age for inclusion

% Automatically determine base directory from JSON file path
json_path = '/mnt/Develop/Ronaldo/dev/Data/NewFolder/XIALPHANET.json';
[dataset_dir, ~, ~] = fileparts(json_path);
% Load and decode dataset JSON
dataset = jsondecode(fileread(json_path));
% Set the location field automatically based on JSON file directory
dataset.Location = dataset_dir;

% Load structural model parameters
parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));
Cortex                  = load("templates/Cortex.mat");




% Positive colormap: white ? orange ? dark orange
n_gray = 95;     % Gray range from 0 to 0.95
n_color = 5;     % Color range from 0.95 to 1
gray = [0.85, 0.85, 0.85];

clip01 = @(x) min(max(x, 0), 1);

% Colors
dark_orange   = clip01([0.60, 0.30, 0.00]);
bright_orange = clip01([1.00, 0.55, 0.00]);

dark_green    = clip01([0.00, 0.35, 0.10]);
bright_green  = clip01([0.10, 0.85, 0.50]);

dark_blue     = clip01([0.15, 0.30, 0.65]);
bright_blue   = clip01([0.30, 0.60, 1.00]);

% Create gray part (flat gray)
gray_part = repmat(gray, n_gray, 1);

% Positive transitions
orange_positive = [gray_part; ...
                   [linspace(gray(1), bright_orange(1), n_color)', ...
                    linspace(gray(2), bright_orange(2), n_color)', ...
                    linspace(gray(3), bright_orange(3), n_color)']];

green_positive = [gray_part; ...
                  [linspace(gray(1), bright_green(1), n_color)', ...
                   linspace(gray(2), bright_green(2), n_color)', ...
                   linspace(gray(3), bright_green(3), n_color)']];

blue_positive = [gray_part; ...
                 [linspace(gray(1), bright_blue(1), n_color)', ...
                  linspace(gray(2), bright_blue(2), n_color)', ...
                  linspace(gray(3), bright_blue(3), n_color)']];

% Negative transitions (also gray to dark color)
orange_negative = [gray_part; ...
                   [linspace(gray(1), dark_orange(1), n_color)', ...
                    linspace(gray(2), dark_orange(2), n_color)', ...
                    linspace(gray(3), dark_orange(3), n_color)']];

green_negative = [gray_part; ...
                  [linspace(gray(1), dark_green(1), n_color)', ...
                   linspace(gray(2), dark_green(2), n_color)', ...
                   linspace(gray(3), dark_green(3), n_color)']];

blue_negative = [gray_part; ...
                 [linspace(gray(1), dark_blue(1), n_color)', ...
                  linspace(gray(2), dark_blue(2), n_color)', ...
                  linspace(gray(3), dark_blue(3), n_color)']];

% Optional preview
figure;
subplot(3,2,1); imagesc(permute(orange_negative, [1 3 2])); title('Orange Negative');
subplot(3,2,2); imagesc(permute(orange_positive, [1 3 2])); title('Orange Positive');

subplot(3,2,3); imagesc(permute(green_negative, [1 3 2])); title('Green Negative');
subplot(3,2,4); imagesc(permute(green_positive, [1 3 2])); title('Green Positive');

subplot(3,2,5); imagesc(permute(blue_negative, [1 3 2]));  title('Blue Negative');
subplot(3,2,6); imagesc(permute(blue_positive, [1 3 2]));  title('Blue Positive');



% Import required helper functions
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*


ages = [];
delays = [];
All_Data = {}; 
index = 1;
N_roi = 360;
for i=1:length(dataset.Participants)
    participant = dataset.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed')) && 0 <= participant_age && participant_age <= 100
       ages = [ages,participant_age];
       All_Data{2,index} =  participant_age;
       Part_Info = jsondecode(fileread(fullfile(dataset.Location,participant.SubID,participant.FileInfo)));
       alpha_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
       a(:,1) = alpha_process.Power;
       a(:,2) = alpha_process.Width;
       a(:,3) = alpha_process.Exponent;
       a(:,4) = alpha_process.PAF;
       xi_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
       e(:,1) = xi_process.Power;
       e(:,2) = xi_process.Width;
       e(:,3) = xi_process.Exponent;
       s2 = 1;
       x = v2x(e,a,s2);
       D = load(fullfile(dataset.Location, participant.SubID, Part_Info.Delay_Matrix));
       delays(index) = 1000 * mean(D.Delay_Matrix(:));
       All_Data{1,index} = x;
       index = index +1;
    end
end

% Initialize storage for Peak Alpha Frequency (PAF), Amplitude of the Alpha, and Amplitude of Xi
% Initialize storage
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

threshold_PAF = 7;
import functions.auxx.Refine_Solution.*

% Temporary storage for thresholds
threshold_Alpha = zeros(1, length(All_Data));
threshold_Xi = zeros(1, length(All_Data));

% First pass: extract features
for j = 1:length(All_Data)
    fprintf('Processing subject %d\n', j);
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Alpha Amplitude
    XiAmp_all(:,j) = e(:,1);          % Xi Amplitude
end

% Compute mean across voxels for each subject
mean_PAF = mean(PAF_all, 1, 'omitnan');
mean_Alpha = mean(AlphaAmp_all, 1, 'omitnan');
mean_Xi = mean(XiAmp_all, 1, 'omitnan');

% Z-score computation
z_PAF = zscore(mean_PAF);
z_Alpha = zscore(mean_Alpha);
z_Xi = zscore(mean_Xi);

% Identify outliers: abs(z) > 3 in any metric
outlier_idx = abs(z_PAF) > 3 | abs(z_Alpha) > 3 | abs(z_Xi) > 3;
valid_idx = ~outlier_idx;
fprintf('Removed %d outlier subjects\n', sum(outlier_idx));

% Filter out outliers
PAF_all = PAF_all(:, valid_idx);
AlphaAmp_all = AlphaAmp_all(:, valid_idx);
XiAmp_all = XiAmp_all(:, valid_idx);
All_Data = All_Data(:, valid_idx);

% Second pass: Apply thresholds
parfor j = 1:length(All_Data)
    [e, a, s2] = x2v(All_Data{1,j});
    Alpha_j = a(:,1);
    Xi_j = e(:,1);
    PAF_j = a(:,4);

    threshold_Alpha(j) = set_threshold_em(Alpha_j);
    threshold_Xi(j) = set_threshold_em(Xi_j);
    if mode == 1 % For probability
        Alpha_j =  (Alpha_j > threshold_Alpha(j));
        Xi_j_bin = Xi_j > threshold_Alpha(j);
        PAF_j_bin = (PAF_j .* (Alpha_j > threshold_Alpha(j))) > threshold_PAF;
    else % For amplitude Distribution
        Alpha_j =  Alpha_j.*(Alpha_j > threshold_Alpha(j));
        Xi_j_bin = Xi_j.*(Xi_j > threshold_Alpha(j));
        PAF_j_bin = PAF_j.*((PAF_j .* (Alpha_j > threshold_Alpha(j))) > threshold_PAF);
    end

    AlphaAmp_all(:,j) = Alpha_j;
    XiAmp_all(:,j) = Xi_j_bin;
    PAF_all(:,j) = PAF_j_bin;
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
delays = delays(valid_idx);
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age

%
delays = delays(:);
N_age = length(delays);

% Compute means across voxels (positive values only)
AlphaAmp_mean = nan(N_age,1);
XiAmp_mean    = nan(N_age,1);
PAF_mean      = nan(N_age,1);

for t = 1:N_age
    a = AlphaAmp_all(:,t); a = a(a > 0);
    x = XiAmp_all(:,t);    x = x(x > 0);
    p = PAF_all(:,t);      p = p(p >= 7 & p<= 13 );  % Only keep PAF > 7

    AlphaAmp_mean(t) = 10*log10(mean(a));  
    XiAmp_mean(t)    = 10*log10(mean(x));
    PAF_mean(t)      = mean(p);
end

% Variable configs: [data, label, color]
variables = {
    AlphaAmp_mean, 'Alpha Amplitude (dB)', orange_negative(end,:);
    XiAmp_mean,    'Xi Amplitude (dB)',    green_negative(end,:);
    PAF_mean,      'PAF (Hz)',             blue_negative(end,:);
};

fontSize = 16;
lineWidth = 2;

for i = 1:3
    y = variables{i,1};
    y_label = variables{i,2};
    c = variables{i,3};

    % --- Remove IQR outliers
    Q1 = quantile(y, 0.25); Q3 = quantile(y, 0.75);
    IQR_val = Q3 - Q1;
    idx = (y >= Q1 - 3*IQR_val) & (y <= Q3 + 3* IQR_val);
    x = delays(idx); 
    y = y(idx);

    % --- Regression
    [x_sorted, idxSort] = sort(x);
    y_sorted = y(idxSort);
    X_quad = [x_sorted];
    [b, stats] = robustfit(X_quad, y_sorted);
    y_fit = b(1) + b(2)*x_sorted;
    pval_b3 = stats.p(2);

    % --- Confidence bands
    X_design = [ones(size(x_sorted)), x_sorted];
    var_fit = sum((X_design * stats.covb) .* X_design, 2);
    se_fit = sqrt(var_fit);
    upper = y_fit + 4*se_fit;
    lower = y_fit - 4*se_fit;

    % === Axes layout setup
    figure('Color','w', 'Units','normalized', 'Position', [0.3 0.3 0.6 0.6])

    ax_main = axes('Position', [0.15 0.15 0.65 0.65]); % main plot
    hold(ax_main, 'on');
    h_fill = fill([x_sorted; flipud(x_sorted)], [upper; flipud(lower)], ...
        c, 'FaceAlpha', 0.25, 'EdgeColor','none', 'Parent', ax_main);
    h_fit = plot(ax_main, x_sorted, y_fit, '-', 'Color', c, 'LineWidth', lineWidth);

    xlabel(ax_main, 'Delays (ms)', 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel(ax_main, y_label, 'FontSize', fontSize, 'FontWeight', 'bold');
    title(ax_main, sprintf('%s vs Delays   (p = %.3g)', y_label, pval_b3), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold', ...
        'Units','normalized', 'Position',[0.5, 1.12, 0]);

    set(ax_main, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlim(ax_main, [min(delays), max(delays)]);  
    grid(ax_main, 'on');

    % === TOP KDE of delays
    ax_top = axes('Position', [0.15 0.76 0.65 0.08]);
    [fx, xgrid] = ksdensity(x_sorted);
    pd_x = fitdist(x_sorted, 'Normal');
    fx_gauss = normpdf(xgrid, pd_x.mu, pd_x.sigma);

    fill(ax_top, xgrid, fx, 'r', ...
        'FaceAlpha', 0.25, 'EdgeColor', 'r', 'LineWidth', 1.5); hold on;
    h_gauss_delay = plot(ax_top, xgrid, fx_gauss, '--', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

    axis(ax_top, 'tight'); set(ax_top, 'XTick', [], 'YTick', []); box off;
    ax_top.XColor = 'none'; ax_top.YColor = 'none';

    % === RIGHT KDE of y
    ax_right = axes('Position', [0.82 0.15 0.12 0.62]);
    [fy, ygrid] = ksdensity(y);
    pd_y = fitdist(y, 'Normal');
    fy_gauss = normpdf(ygrid, pd_y.mu, pd_y.sigma);

    fill(ax_right, fy, ygrid, c, ...
        'FaceAlpha', 0.25, 'EdgeColor', c, 'LineWidth', 1.5); hold on;
    h_gauss_y = plot(ax_right, fy_gauss, ygrid, '--', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

    axis(ax_right, 'tight'); set(ax_right, 'XTick', [], 'YTick', []); box off;
    ax_right.XColor = 'none'; ax_right.YColor = 'none';

    % === Legend in main axis
    legend(ax_main, [h_fit, h_fill, h_gauss_y], ...
        {'Regression Fit', 'Estimator Uncertainty', 'Gaussian Fit (Y)'}, ...
        'Location', 'southwest', 'FontSize', 12);
end

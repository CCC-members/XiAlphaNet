clc;
clear all;

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

% Path to the JSON file with model result metadata. Modify this directions
% manually acording to the location of the downloaded data 
json_path = '/mnt/Develop/Ronaldo/program_working/xialphanet_newresults22/XIALPHANET.json';
dir_data = '/mnt/Develop/Ronaldo/dev/Data/norms';

% Automatically determine base directory from JSON file path
[dataset_dir, ~, ~] = fileparts(json_path);
% Load and decode dataset JSON
dataset = jsondecode(fileread(json_path));
% Set the location field automatically based on JSON file directory
dataset.Location = dataset_dir;

% Load structural model parameters
parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));
Cortex                  = load("templates/Cortex.mat");




%% Positive colormap: white ? orange ? dark orange
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



%%
% Import required helper functions
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*


ages = [];
All_Data = {}; 
index = 1;
HCP_MMP1 = Cortex.Atlas(Cortex.iAtlas);
N_voxel = length(Cortex.Curvature);
N_roi = 360;
for i=1:length(dataset.Participants)
    participant = dataset.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed')) && age_min <= participant_age && participant_age <= age_max
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
ages(isnan(ages)) = mean(ages(~isnan(ages)));  % Replace NaNs with mean age


%% === COMPONENT CONFIGURATION ===
components     = {'PAF', 'AlphaAmp', 'XiAmp'};
data_all       = {PAF_all, AlphaAmp_all, XiAmp_all};
colormaps_all  = {
    orange_positive, orange_negative; 
    blue_positive,   blue_negative; 
    green_positive,  green_negative
};
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*

% === PROCESS EACH COMPONENT (VOXEL-WISE) ===
for comp_idx = 1:3
    fprintf('\nProcessing %s (voxel-wise)\n', components{comp_idx});
    current_data = data_all{comp_idx};
    cmap_pos     = colormaps_all{comp_idx, 1};
    cmap_neg     = colormaps_all{comp_idx, 2};

    % === Fit ZIG model for each voxel ===
    result_all = cell(N_voxel, 1);
    for v = 1:N_voxel
        v
        Y_v = current_data(v, :);       % [1 x N] vector for voxel v
        result_all{v} = fitZIG_random(Y_v, ages);
    end

    % === GAMMA > 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.gamma(2) > 0
            J(v) = 1 - result_all{v}.p_gamma(2);
        end
    end
    guide.Visualization.esi_plot_single2;
    colormap(cmap_pos); caxis([0 1]); title([components{comp_idx}, ' \gamma > 0']);

    % === GAMMA < 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.gamma(2) < 0
            J(v) = 1 - result_all{v}.p_gamma(2);
        end
    end
    guide.Visualization.esi_plot_single2;
    colormap(cmap_neg); caxis([0 1]); title([components{comp_idx}, ' \gamma < 0']);

    % === BETA > 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.beta(3) > 0
            J(v) = 1 - result_all{v}.p_beta(3);
        end
    end
    guide.Visualization.esi_plot_single2;
    colormap(cmap_pos); caxis([0 1]); title([components{comp_idx}, ' \beta > 0']);

    % === BETA < 0 ===
    J = zeros(N_voxel, 1);
    for v = 1:N_voxel
        if result_all{v}.beta(3) < 0
            J(v) = 1 - result_all{v}.p_beta(3);
        end
    end
    guide.Visualization.esi_plot_single2;
    colormap(cmap_neg); caxis([0 1]); title([components{comp_idx}, ' \beta < 0']);
end

%%
ages = ages(:);
N_age = length(ages);

% === Compute means across voxels (positive values only) ===
AlphaAmp_mean = nan(N_age,1);
XiAmp_mean    = nan(N_age,1);
PAF_mean      = nan(N_age,1);

for t = 1:N_age
    a = AlphaAmp_all(:,t); a = a(a > 0);
    x = XiAmp_all(:,t);    x = x(x > 0);
    p = PAF_all(:,t);      p = p(p > 0);

    AlphaAmp_mean(t) = 10*log10(mean(a));  
    XiAmp_mean(t)    = 10*log10(mean(x));
    PAF_mean(t)      = mean(p);
end

% === Variable configs: [data, label, color] ===
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

    % --- Remove IQR outliers ---
    Q1 = quantile(y, 0.25); Q3 = quantile(y, 0.75);
    IQR_val = Q3 - Q1;
    idx = (y >= Q1 - 1.5*IQR_val) & (y <= Q3 + 1.5*IQR_val);
    x = ages(idx); 
    y = y(idx);

    % --- Regression ---
    [x_sorted, idxSort] = sort(x);
    y_sorted = y(idxSort);
    X_quad = [x_sorted, x_sorted.^2];
    [b, stats] = robustfit(X_quad, y_sorted);
    y_fit = b(1) + b(2)*x_sorted + b(3)*x_sorted.^2;
    pval_b3 = stats.p(3);

    % --- Classical R² and f² ---
    y_mean = mean(y_sorted);
    SSR = sum((y_sorted - y_fit).^2);
    SST = sum((y_sorted - y_mean).^2);
    R2 = 1 - SSR/SST;
    f2 = R2 / (1 - R2);

    % --- Confidence and SE bands ---
    X_design = [ones(size(x_sorted)), x_sorted, x_sorted.^2];
    var_fit = sum((X_design * stats.covb) .* X_design, 2);
    se_fit = sqrt(var_fit);
    upper_CI = y_fit + 2*se_fit;
    lower_CI = y_fit - 2*se_fit;
    upper_SE = y_fit + se_fit;
    lower_SE = y_fit - se_fit;

    % === FIGURE ===
    fig = figure('Color','w', 'Units','normalized', 'Position', [0.3 0.3 0.6 0.6]);

    % --- Main plot ---
    ax_main = axes('Parent', fig, 'Position', [0.15 0.15 0.65 0.65]);
    hold(ax_main, 'on');

    % 95% CI (lighter)
    fill([x_sorted; flipud(x_sorted)], [upper_CI; flipud(lower_CI)], ...
        c, 'FaceAlpha', 0.15, 'EdgeColor','none', 'Parent', ax_main);
    % ±1 SE (darker)
    fill([x_sorted; flipud(x_sorted)], [upper_SE; flipud(lower_SE)], ...
        c, 'FaceAlpha', 0.25, 'EdgeColor','none', 'Parent', ax_main);

    plot(ax_main, x_sorted, y_fit, '-', 'Color', c, 'LineWidth', lineWidth);

    % === Format p-value string ===
    if pval_b3 < 0.001
        p_str = 'p < 0.001';
    else
        p_str = sprintf('p = %.3f', pval_b3);
    end

    % === Labels ===
    xlabel(ax_main, sprintf('Age (Years, R² = %.3f, f² = %.3f, %s)', R2, f2, p_str), ...
        'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel(ax_main, y_label, 'FontSize', fontSize, 'FontWeight', 'bold');
    set(ax_main, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlim(ax_main, [5 95]); grid(ax_main, 'on');

    % === Top KDE ===
    ax_top = axes('Parent', fig, 'Position', [0.15 0.81 0.65 0.12]);
    [fx, xgrid] = ksdensity(x_sorted);
    pd_x = fitdist(x_sorted, 'Normal');
    fx_gauss = normpdf(xgrid, pd_x.mu, pd_x.sigma);

    fill(ax_top, xgrid, fx, [0.5 0.5 0.5], ...
        'FaceAlpha', 0.25, 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
    hold on;
    plot(ax_top, xgrid, fx_gauss, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    axis(ax_top, 'tight'); set(ax_top, 'XTick', [], 'YTick', []); box off;
    xlim(ax_top, [5 95]);
    ax_top.XColor = 'none'; ax_top.YColor = 'none';

    % === Right KDE ===
    ax_right = axes('Parent', fig, 'Position', [0.82 0.15 0.12 0.65]);
    [fy, ygrid] = ksdensity(y);
    pd_y = fitdist(y, 'Normal');
    fy_gauss = normpdf(ygrid, pd_y.mu, pd_y.sigma);

    fill(ax_right, fy, ygrid, c, ...
        'FaceAlpha', 0.25, 'EdgeColor', c, 'LineWidth', 1.5);
    hold on;
    plot(ax_right, fy_gauss, ygrid, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    axis(ax_right, 'tight'); set(ax_right, 'XTick', [], 'YTick', []); box off;
    ax_right.XColor = 'none'; ax_right.YColor = 'none';

    % === Top title via annotation (always visible and centered) ===
    annotation(fig, 'textbox', [0.15, 0.92, 0.65, 0.05], ...
        'String', sprintf('%s vs Age', y_label), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', fontSize+2, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

%%

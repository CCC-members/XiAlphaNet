%Delayalpha
% --------------------------- Setup ----------------------------------

% Clear workspace and command window for a fresh start
clear; clc;


%Directory containing .mat files

% Path to the JSON file with model result metadata
json_path = '/mnt/Store/Ronaldo/dev/Data/NewFolder/XIALPHANET.json';
% Automatically determine base directory from JSON file path
[dataset_dir, ~, ~] = fileparts(json_path);
% Load and decode dataset JSON
dataset = jsondecode(fileread(json_path));

% Set the location field automatically based on JSON file directory
dataset.Location = dataset_dir;
import templates.*
import functions.auxx.ModelVectorization.* 
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*
load("templates/mylin_data.mat")
% Initialize arrays for storing delays and ages
delays = [];
ages = [];
alpha_powers= [];
xi_powers= [];
alpha_pfs = [];


%--------------------------- Data Extraction -----------------------
index = 1;
parfor i=1:length(dataset.Participants)
    i
    participant = dataset.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed'))
        ages(i) = participant_age;
        Part_Info = jsondecode(fileread(fullfile(dataset.Location,participant.SubID,participant.FileInfo)));
        Mod_Weights = load(fullfile(dataset.Location,participant.SubID,Part_Info.Mod_Weights));
        Delay_Matrix = load(fullfile(dataset.Location,participant.SubID,Part_Info.Delay_Matrix));

        Alpha_estimate = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
        Xi_estimate = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
        %
        threshold_alpha =  set_threshold_em(Alpha_estimate.Power);
        pos_alpha = (Alpha_estimate.Power>threshold_alpha);
        alpha_powers(i) =real( mean(Alpha_estimate.Power(pos_alpha)));
        alpha_pfs(i) = real(mean(Alpha_estimate.PAF(pos_alpha)));
        threshold_xi = set_threshold_em(Xi_estimate.Power);
        pos_xi = (Xi_estimate.Power>threshold_xi);
        xi_powers(i) = real(mean(Xi_estimate.Power(pos_xi)));
        delays(i) =  mean(Delay_Matrix.Delay_Matrix(:));
        index = index +1;
    end
end



%%
delays0 = delays(:);
ages0 = ages(:);
alpha_powers0 = alpha_powers(:);
alpha_pfs0 = alpha_pfs(:);
xi_powers0 = xi_powers(:);
%%
delays = delays0(:);
ages = ages0(:);
alpha_powers = alpha_powers0(:);
alpha_pfs = alpha_pfs0(:);
xi_powers = xi_powers0(:);


%%
%% --------------------------- Alpha Power vs Delays ---------------------------
% 1) Log-transform and sort
AP = 10*log10(1e-6 + alpha_powers);          % log(Alpha Power)
%AP =  (AP(:)-min(AP(:)))./(max(AP(:))-min(AP(:)));
Del = 1000*(delays0);               % log(tau^-2)
%Del =  (Del(:)-min(Del(:)))./(max(Del(:))-min(Del(:)));
% Example IQR-based outlier removal on PAF:
Q1 = quantile(AP, 0.25);
Q3 = quantile(AP, 0.75);
IQR_val = Q3 - Q1;

lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (AP < lowerBound) | (AP > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

AP    = AP(~idxOut);
Del = Del(~idxOut);


[logDel_sorted, idxSort] = sort(Del);
logAP_sorted = AP(idxSort);

% 2) Robust Linear Fit
X_lin = logDel_sorted;
[b_lin, stats_lin] = robustfit(X_lin, logAP_sorted);

% Evaluate linear fit on the sorted data
y_lin_fit = b_lin(1) + b_lin(2) * logDel_sorted;

% 1-StdErr confidence bands for linear fit
X_lin_fit = [ones(size(logDel_sorted)), logDel_sorted];
var_lin_fit = sum((X_lin_fit * stats_lin.covb) .* X_lin_fit, 2);
se_lin_fit = sqrt(var_lin_fit);
upper_lin = y_lin_fit + se_lin_fit;
lower_lin = y_lin_fit - se_lin_fit;

% 3) Robust Quadratic Fit
X_quad = [logDel_sorted, logDel_sorted.^2];
[b_quad, stats_quad] = robustfit(X_quad, logAP_sorted);

y_quad_fit = b_quad(1) + b_quad(2) * logDel_sorted + b_quad(3) * logDel_sorted.^2;

% 1-StdErr confidence bands for quadratic fit
X_quad_fit = [ones(size(logDel_sorted)), logDel_sorted, logDel_sorted.^2];
var_quad_fit = sum((X_quad_fit * stats_quad.covb) .* X_quad_fit, 2);
se_quad_fit = sqrt(var_quad_fit);
upper_quad = y_quad_fit + se_quad_fit;
lower_quad = y_quad_fit - se_quad_fit;

% 4) Model comparison (RMSE and BIC)
% Linear model residuals and RMSE
residuals_lin = logAP_sorted - y_lin_fit;
rmse_lin = sqrt(mean(residuals_lin.^2));

% Quadratic model residuals and RMSE
residuals_quad = logAP_sorted - y_quad_fit;
rmse_quad = sqrt(mean(residuals_quad.^2));

% BIC for Linear Model
n = length(logAP_sorted);
k_lin = length(b_lin); % Number of parameters in linear model
bic_lin = n*log(mean(residuals_lin.^2)) + k_lin*log(n);

% BIC for Quadratic Model
k_quad = length(b_quad); % Number of parameters in quadratic model
bic_quad = n*log(mean(residuals_quad.^2)) + k_quad*log(n);

% 5) Plot only the regression lines (no scatter)
figure; hold on;

% Linear fit (red line)
plot(logDel_sorted, y_lin_fit, 'r-', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_lin; flipud(lower_lin)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Quadratic fit (blue dashed)
plot(logDel_sorted, y_quad_fit, 'b--', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_quad; flipud(lower_quad)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Axis labels are logs of the variables
xlabel('Delays (ms)', 'FontWeight', 'bold');
ylabel('Alpha Power Amplitud (dB)', 'FontWeight', 'bold');
title('Average Source Alpha Power ~ \tau  (Robust Linear & Quadratic)');
xlim([4 18])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;
%% --------------------------- Peak Alpha Frequency vs Delays ---------------------------
% 1) Log-transform and sort
AP = alpha_pfs; % log(Alpha Power)
%AP =  (AP(:)-min(AP(:)))./(max(AP(:))-min(AP(:)));
Del = (1000*delays0);               % log(tau^-2)
%Del =  (Del(:)-min(Del(:)))./(max(Del(:))-min(Del(:)));

% Example IQR-based outlier removal on PAF:
Q1 = quantile(AP, 0.25);
Q3 = quantile(AP, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (AP < lowerBound) | (AP > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

AP    = AP(~idxOut);
Del = Del(~idxOut);




[logDel_sorted, idxSort] = sort(Del);
logAP_sorted = AP(idxSort);

% 2) Robust Linear Fit
X_lin = logDel_sorted;
[b_lin, stats_lin] = robustfit(X_lin, logAP_sorted);

% Evaluate linear fit on the sorted data
y_lin_fit = b_lin(1) + b_lin(2) * logDel_sorted;

% 1-StdErr confidence bands for linear fit
X_lin_fit = [ones(size(logDel_sorted)), logDel_sorted];
var_lin_fit = sum((X_lin_fit * stats_lin.covb) .* X_lin_fit, 2);
se_lin_fit = sqrt(var_lin_fit);
upper_lin = y_lin_fit + se_lin_fit;
lower_lin = y_lin_fit - se_lin_fit;

% 3) Robust Quadratic Fit
X_quad = [logDel_sorted, logDel_sorted.^2];
[b_quad, stats_quad] = robustfit(X_quad, logAP_sorted);

y_quad_fit = b_quad(1) + b_quad(2) * logDel_sorted + b_quad(3) * logDel_sorted.^2;

% 1-StdErr confidence bands for quadratic fit
X_quad_fit = [ones(size(logDel_sorted)), logDel_sorted, logDel_sorted.^2];
var_quad_fit = sum((X_quad_fit * stats_quad.covb) .* X_quad_fit, 2);
se_quad_fit = sqrt(var_quad_fit);
upper_quad = y_quad_fit + se_quad_fit;
lower_quad = y_quad_fit - se_quad_fit;

% 4) Model comparison (RMSE and BIC)
% Linear model residuals and RMSE
residuals_lin = logAP_sorted - y_lin_fit;
rmse_lin = sqrt(mean(residuals_lin.^2));

% Quadratic model residuals and RMSE
residuals_quad = logAP_sorted - y_quad_fit;
rmse_quad = sqrt(mean(residuals_quad.^2));

% BIC for Linear Model
n = length(logAP_sorted);
k_lin = length(b_lin); % Number of parameters in linear model
bic_lin = n*log(mean(residuals_lin.^2)) + k_lin*log(n);

% BIC for Quadratic Model
k_quad = length(b_quad); % Number of parameters in quadratic model
bic_quad = n*log(mean(residuals_quad.^2)) + k_quad*log(n);

% 5) Plot only the regression lines (no scatter)
figure; hold on;

% Linear fit (red line)
plot(logDel_sorted, y_lin_fit, 'r-', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_lin; flipud(lower_lin)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Quadratic fit (blue dashed)
plot(logDel_sorted, y_quad_fit, 'b--', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_quad; flipud(lower_quad)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Axis labels are logs of the variables
xlabel('Delays (ms)', 'FontWeight', 'bold');
ylabel('PAF (Hz)', 'FontWeight', 'bold');
title('Average Source Peak Alpha Frequency (PAF) ~ \tau  (Robust Linear & Quadratic)');
%xlim([4 18])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

%% --------------------------- Xi Power vs Delays ---------------------------
XIP = 10*log10(10^(-6)+xi_powers); 
%XIP =  (XIP(:)-min(XIP(:)))./(max(XIP(:))-min(XIP(:)));
Del = 1000*(delays0);               % log(tau^-2)
%Del =  (Del(:)-min(Del(:)))./(max(Del(:))-min(Del(:)));

% Example IQR-based outlier removal on PAF:
Q1 = quantile(XIP, 0.25);
Q3 = quantile(XIP, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (XIP < lowerBound) | (XIP > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

XIP    = XIP(~idxOut);
Del = Del(~idxOut);


[logDel_sorted, idxSort] = sort(Del);
logXIP_sorted = XIP(idxSort);

% Robust Linear
X_lin = logDel_sorted;
[b_lin, stats_lin] = robustfit(X_lin, logXIP_sorted);

y_lin_fit = b_lin(1) + b_lin(2) * logDel_sorted;

% 1-StdErr confidence bands for linear fit
X_lin_fit = [ones(size(logDel_sorted)), logDel_sorted];
var_lin_fit = sum((X_lin_fit * stats_lin.covb) .* X_lin_fit, 2);
se_lin_fit = sqrt(var_lin_fit);
upper_lin = y_lin_fit + se_lin_fit;
lower_lin = y_lin_fit - se_lin_fit;

% Robust Quadratic
X_quad = [logDel_sorted, logDel_sorted.^2];
[b_quad, stats_quad] = robustfit(X_quad, logXIP_sorted);

y_quad_fit = b_quad(1) + b_quad(2) * logDel_sorted + b_quad(3) * logDel_sorted.^2;

% 1-StdErr confidence bands for quadratic fit
X_quad_fit = [ones(size(logDel_sorted)), logDel_sorted, logDel_sorted.^2];
var_quad_fit = sum((X_quad_fit * stats_quad.covb) .* X_quad_fit, 2);
se_quad_fit = sqrt(var_quad_fit);
upper_quad = y_quad_fit + se_quad_fit;
lower_quad = y_quad_fit - se_quad_fit;

% Model comparison (RMSE and BIC)
% Linear model residuals and RMSE
residuals_lin = logXIP_sorted - y_lin_fit;
rmse_lin = sqrt(mean(residuals_lin.^2));

% Quadratic model residuals and RMSE
residuals_quad = logXIP_sorted - y_quad_fit;
rmse_quad = sqrt(mean(residuals_quad.^2));

% BIC for Linear Model
n = length(logXIP_sorted);
k_lin = length(b_lin); % Number of parameters in linear model
bic_lin = n*log(mean(residuals_lin.^2)) + k_lin*log(n);

% BIC for Quadratic Model
k_quad = length(b_quad); % Number of parameters in quadratic model
bic_quad = n*log(mean(residuals_quad.^2)) + k_quad*log(n);

% Plot
figure; hold on;

plot(logDel_sorted, y_lin_fit, 'r-', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_lin; flipud(lower_lin)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(logDel_sorted, y_quad_fit, 'b--', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_quad; flipud(lower_quad)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Delay (ms)', 'FontWeight', 'bold');
ylabel('Xi Power Amplitud (dB)', 'FontWeight', 'bold');
title('Average Source Xi Power ~ \tau (Robust Linear & Quadratic)');
xlim([4 18])
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

%% --------------------------- Alpha Power vs Age ---------------------------
% 1) Log-transform and sort
AP = 10*log10(1e-6 + alpha_powers);          % log(Alpha Power)
%AP =  (AP(:)-min(AP(:)))./(max(AP(:))-min(AP(:)));
Del = ages;               % log(tau^-2)
% Example IQR-based outlier removal on PAF:
Q1 = quantile(AP, 0.25);
Q3 = quantile(AP, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (AP < lowerBound) | (AP > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

AP    = AP(~idxOut);
Del = Del(~idxOut);


[logDel_sorted, idxSort] = sort(Del);
logAP_sorted = AP(idxSort);

% 2) Robust Linear Fit
X_lin = logDel_sorted;
[b_lin, stats_lin] = robustfit(X_lin, logAP_sorted);

% Evaluate linear fit on the sorted data
y_lin_fit = b_lin(1) + b_lin(2) * logDel_sorted;

% 1-StdErr confidence bands for linear fit
X_lin_fit = [ones(size(logDel_sorted)), logDel_sorted];
var_lin_fit = sum((X_lin_fit * stats_lin.covb) .* X_lin_fit, 2);
se_lin_fit = sqrt(var_lin_fit);
upper_lin = y_lin_fit + 2*se_lin_fit;
lower_lin = y_lin_fit - 2*se_lin_fit;

% 3) Robust Quadratic Fit
X_quad = [logDel_sorted, logDel_sorted.^2];
[b_quad, stats_quad] = robustfit(X_quad, logAP_sorted);

y_quad_fit = b_quad(1) + b_quad(2) * logDel_sorted + b_quad(3) * logDel_sorted.^2;

% 1-StdErr confidence bands for quadratic fit
X_quad_fit = [ones(size(logDel_sorted)), logDel_sorted, logDel_sorted.^2];
var_quad_fit = sum((X_quad_fit * stats_quad.covb) .* X_quad_fit, 2);
se_quad_fit = sqrt(var_quad_fit);
upper_quad = y_quad_fit + 2*se_quad_fit;
lower_quad = y_quad_fit - 2*se_quad_fit;

% 4) Model comparison (RMSE and BIC)
% Linear model residuals and RMSE
residuals_lin = logAP_sorted - y_lin_fit;
rmse_lin = sqrt(mean(residuals_lin.^2));

% Quadratic model residuals and RMSE
residuals_quad = logAP_sorted - y_quad_fit;
rmse_quad = sqrt(mean(residuals_quad.^2));

% BIC for Linear Model
n = length(logAP_sorted);
k_lin = length(b_lin); % Number of parameters in linear model
bic_lin = n*log(mean(residuals_lin.^2)) + k_lin*log(n);

% BIC for Quadratic Model
k_quad = length(b_quad); % Number of parameters in quadratic model
bic_quad = n*log(mean(residuals_quad.^2)) + k_quad*log(n);

% 5) Plot only the regression lines (no scatter)
figure; hold on;

% Linear fit (red line)
plot(logDel_sorted, y_lin_fit, 'r-', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_lin; flipud(lower_lin)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Quadratic fit (blue dashed)
plot(logDel_sorted, y_quad_fit, 'b--', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_quad; flipud(lower_quad)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Axis labels are logs of the variables
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('Alpha Power Amplitud (dB)', 'FontWeight', 'bold');
title('Average Source Alpha Power ~ age  (Robust Linear & Quadratic)');
xlim([5 95])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;
% --------------------------- Peak Alpha Frequency vs Age ---------------------------
% 1) Log-transform and sort
AP = alpha_pfs; % log(Alpha Power)
%AP =  (AP(:)-min(AP(:)))./(max(AP(:))-min(AP(:)));
Del = ages;              % log(tau^-2)

% Example IQR-based outlier removal on PAF:
Q1 = quantile(AP, 0.25);
Q3 = quantile(AP, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (AP < lowerBound) | (AP > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

AP    = AP(~idxOut);
Del = Del(~idxOut);




[logDel_sorted, idxSort] = sort(Del);
logAP_sorted = AP(idxSort);

% 2) Robust Linear Fit
X_lin = logDel_sorted;
[b_lin, stats_lin] = robustfit(X_lin, logAP_sorted);

% Evaluate linear fit on the sorted data
y_lin_fit = b_lin(1) + b_lin(2) * logDel_sorted;

% 1-StdErr confidence bands for linear fit
X_lin_fit = [ones(size(logDel_sorted)), logDel_sorted];
var_lin_fit = sum((X_lin_fit * stats_lin.covb) .* X_lin_fit, 2);
se_lin_fit = sqrt(var_lin_fit);
upper_lin = y_lin_fit + 2*se_lin_fit;
lower_lin = y_lin_fit - 2*se_lin_fit;

% 3) Robust Quadratic Fit
X_quad = [logDel_sorted, logDel_sorted.^2];
[b_quad, stats_quad] = robustfit(X_quad, logAP_sorted);

y_quad_fit = b_quad(1) + b_quad(2) * logDel_sorted + b_quad(3) * logDel_sorted.^2;

% 1-StdErr confidence bands for quadratic fit
X_quad_fit = [ones(size(logDel_sorted)), logDel_sorted, logDel_sorted.^2];
var_quad_fit = sum((X_quad_fit * stats_quad.covb) .* X_quad_fit, 2);
se_quad_fit = sqrt(var_quad_fit);
upper_quad = y_quad_fit + 2*se_quad_fit;
lower_quad = y_quad_fit - 2*se_quad_fit;

% 4) Model comparison (RMSE and BIC)
% Linear model residuals and RMSE
residuals_lin = logAP_sorted - y_lin_fit;
rmse_lin = sqrt(mean(residuals_lin.^2));

% Quadratic model residuals and RMSE
residuals_quad = logAP_sorted - y_quad_fit;
rmse_quad = sqrt(mean(residuals_quad.^2));

% BIC for Linear Model
n = length(logAP_sorted);
k_lin = length(b_lin); % Number of parameters in linear model
bic_lin = n*log(mean(residuals_lin.^2)) + k_lin*log(n);

% BIC for Quadratic Model
k_quad = length(b_quad); % Number of parameters in quadratic model
bic_quad = n*log(mean(residuals_quad.^2)) + k_quad*log(n);

% 5) Plot only the regression lines (no scatter)
figure; hold on;

% Linear fit (red line)
plot(logDel_sorted, y_lin_fit, 'r-', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_lin; flipud(lower_lin)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Quadratic fit (blue dashed)
plot(logDel_sorted, y_quad_fit, 'b--', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_quad; flipud(lower_quad)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Axis labels are logs of the variables
xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('PAF (Hz)', 'FontWeight', 'bold');
title('Average Source Peak Alpha Frequency (PAF) ~ age  (Robust Linear & Quadratic)');
xlim([5 95])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

% --------------------------- Xi Power vs Age ---------------------------
XIP = 10*log10(10^(-6)+xi_powers); 
%XIP =  (XIP(:)-min(XIP(:)))./(max(XIP(:))-min(XIP(:)));
Del = ages;               % log(tau^-2)


% Example IQR-based outlier removal on PAF:
Q1 = quantile(XIP, 0.25);
Q3 = quantile(XIP, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (XIP < lowerBound) | (XIP > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

XIP    = XIP(~idxOut);
Del = Del(~idxOut);



[logDel_sorted, idxSort] = sort(Del);
logXIP_sorted = XIP(idxSort);

% Robust Linear
X_lin = logDel_sorted;
[b_lin, stats_lin] = robustfit(X_lin, logXIP_sorted);

y_lin_fit = b_lin(1) + b_lin(2) * logDel_sorted;

% 1-StdErr confidence bands for linear fit
X_lin_fit = [ones(size(logDel_sorted)), logDel_sorted];
var_lin_fit = sum((X_lin_fit * stats_lin.covb) .* X_lin_fit, 2);
se_lin_fit = sqrt(var_lin_fit);
upper_lin = y_lin_fit + 2*se_lin_fit;
lower_lin = y_lin_fit - 2*se_lin_fit;

% Robust Quadratic
X_quad = [logDel_sorted, logDel_sorted.^2];
[b_quad, stats_quad] = robustfit(X_quad, logXIP_sorted);

y_quad_fit = b_quad(1) + b_quad(2) * logDel_sorted + b_quad(3) * logDel_sorted.^2;

% 1-StdErr confidence bands for quadratic fit
X_quad_fit = [ones(size(logDel_sorted)), logDel_sorted, logDel_sorted.^2];
var_quad_fit = sum((X_quad_fit * stats_quad.covb) .* X_quad_fit, 2);
se_quad_fit = sqrt(var_quad_fit);
upper_quad = y_quad_fit + 2*se_quad_fit;
lower_quad = y_quad_fit - 2*se_quad_fit;

% Model comparison (RMSE and BIC)
% Linear model residuals and RMSE
residuals_lin = logXIP_sorted - y_lin_fit;
rmse_lin = sqrt(mean(residuals_lin.^2));

% Quadratic model residuals and RMSE
residuals_quad = logXIP_sorted - y_quad_fit;
rmse_quad = sqrt(mean(residuals_quad.^2));

% BIC for Linear Model
n = length(logXIP_sorted);
k_lin = length(b_lin); % Number of parameters in linear model
bic_lin = n*log(mean(residuals_lin.^2)) + k_lin*log(n);

% BIC for Quadratic Model
k_quad = length(b_quad); % Number of parameters in quadratic model
bic_quad = n*log(mean(residuals_quad.^2)) + k_quad*log(n);

% Plot
figure; hold on;

plot(logDel_sorted, y_lin_fit, 'r-', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_lin; flipud(lower_lin)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(logDel_sorted, y_quad_fit, 'b--', 'LineWidth', 1.5);
fill([logDel_sorted; flipud(logDel_sorted)], ...
     [upper_quad; flipud(lower_quad)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Age (Years)', 'FontWeight', 'bold');
ylabel('Xi Power Amplitud (dB)', 'FontWeight', 'bold');
title('Average Source Xi Power ~ age (Robust Linear & Quadratic)');
xlim([5 95])
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

% Outlier Elimination and Quadratic Regression Script
% This script reads .mat files from a specified folder, extracts conduction delay and age,
% handles NaN values, eliminates outliers, performs quadratic regression, and plots the results.

% --------------------------- Setup ----------------------------------

% Clear workspace and command window for a fresh start
clear; clc;

%Directory containing .mat files
dataset = jsondecode(fileread('/mnt/Store/Ronaldo/dev/Data/NewFolder/XIALPHANET.json'));
dataset.Location = '/mnt/Store/Ronaldo/dev/Data/NewFolder';
import templates.*
load("templates/mylin_data.mat")
%Initialize arrays for storing delays and ages
delays = [];
ages = [];
%--------------------------- Data Extraction -----------------------
index = 1;
for i=1:length(dataset.Participants)
    participant = dataset.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed'))
        ages(index) = participant_age;
        Part_Info = jsondecode(fileread(fullfile(dataset.Location,participant.SubID,participant.FileInfo)));
        Delay_Matrix = load(fullfile(dataset.Location,participant.SubID,Part_Info.Delay_Matrix));
        delays(index) =  1000*mean(Delay_Matrix.Delay_Matrix(:));
        index = index +1;
    end
end


% Convert to column vectors
delays0=delays(:);
delays = delays(:);
ages = ages(:);

%%
% Remove any remaining NaN values (precautionary step)
validIdx = ~isnan(delays) & ~isnan(ages);
delays = delays(validIdx);
ages = ages(validIdx);

% Check if there are enough valid data points
if length(delays) < 10
    error('Not enough valid data points for regression. Found %d valid points.', length(delays));
end

% --------------------------- Outlier Detection and Removal -------

% Choose an outlier detection method
% Options:
% 1. Z-Score Method
% 2. Interquartile Range (IQR) Method

% Select the method by setting methodVariable to 'zscore' or 'iqr'
method = 'zscore'; % Change to 'iqr' if preferred

switch lower(method)
    case 'zscore'
        % ------------------ Z-Score Outlier Removal -----------------
        % Calculate Z-scores for delays
        zScores = abs(zscore(delays));

        % Define a Z-score threshold (commonly 3)
        zThreshold = 3;

        % Identify non-outlier indices
        nonOutliers = zScores < zThreshold;

        % Display the number of outliers detected
        numOutliers = sum(~nonOutliers);
        fprintf('Z-Score Method: Detected and removing %d outlier(s).\n', numOutliers);

    case 'iqr'
        % ------------------ IQR Outlier Removal ---------------------
        % Calculate the first and third quartiles
        Q1 = prctile(delays, 25);
        Q3 = prctile(delays, 75);
        IQR = Q3 - Q1;

        % Define outlier bounds
        lowerBound = Q1 - 2.5 * IQR;
        upperBound = Q3 + 2.5 * IQR;

        % Identify non-outlier indices
        nonOutliers = (delays >= lowerBound) & (delays <= upperBound);

        % Display the number of outliers detected
        numOutliers = sum(~nonOutliers);
        fprintf('IQR Method: Detected and removing %d outlier(s).\n', numOutliers);

    otherwise
        error('Invalid outlier detection method selected. Choose "zscore" or "iqr".');
end

% Apply outlier removal
cleanDelays = delays(nonOutliers);
cleanAges = ages(nonOutliers);


% Display the number of data points before and after outlier removal
fprintf('Data Points Before Outlier Removal: %d\n', length(delays));
fprintf('Data Points After Outlier Removal: %d\n', length(cleanDelays));

% Ensure there are enough data points after outlier removal
if length(cleanDelays) < 10
    error('Not enough data points after outlier removal for regression. Found %d points.', length(cleanDelays));
end

% --------------------------- Quadratic Regression ----------------

% Sort data by age for proper plotting
[ages, sortIdx] = sort(cleanAges);
delays = cleanDelays(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(ages)), ages, ages.^2];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), delays); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
ages_fit = linspace(min(ages), max(ages), 200)'; % Same range as LOESS fit
robust_fit_interp = b(1) + b(2) * ages_fit + b(3) * ages_fit.^2;
% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(ages_fit)), ages_fit, ages_fit.^2];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  se_robust_fit;
lower_bound_robust = robust_fit_interp -  se_robust_fit;

% Plot the robust quadratic fit in red
plot(ages_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([ages_fit; flipud(ages_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Age');
ylabel('Conduction Delay');
title('Lifespan of Conductions Delays (\tau)');
legend('Estimated Conductions Delays','Estimator Uncertainty');
grid on;
hold off;
stats.p
% --------------------------- Informative Message -------------------

% Display a message with the number of data points used
fprintf('Quadratic regression performed on %d valid data points after outlier removal.\n', length(cleanDelays));
%%
% Define the ages to fit
% Construct design matrix for quadratic robust fit
invsquare_delays = log(1./(cleanDelays.^2));
ages = cleanAges;
pos = ages< 83;
ages = ages(pos);
invsquare_delays = invsquare_delays(pos);

X = [ones(size(ages)), ages, ages.^2];
%
 %Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), invsquare_delays); % Exclude intercept from X in robustfit

% Visualize robust regression results

% Evaluate the robust quadratic fit on the same range as ages_fit
ages_fit = linspace(0.3, 83, 100)'; % Same range as LOESS fit
robust_fit_interp = b(1) + b(2) * ages_fit + b(3) * ages_fit.^2;
robust_fit_interp_exp = exp(b(1) + b(2) * ages_fit + b(3) * ages_fit.^2);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(ages_fit)), ages_fit, ages_fit.^2];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = exp(robust_fit_interp +  se_robust_fit);
lower_bound_robust = exp(robust_fit_interp -  se_robust_fit);

inv_delays = (robust_fit_interp_exp - min(lower_bound_robust)) / (max(upper_bound_robust) - min(lower_bound_robust));
upper_delays = (upper_bound_robust - min(lower_bound_robust)) / (max(upper_bound_robust) - min(lower_bound_robust));
lower_delays = (lower_bound_robust - min(lower_bound_robust)) / (max(upper_bound_robust) - min(lower_bound_robust));

% Plot the results
figure;
hold on;


% Plot the robust quadratic fit in red
plot(ages_fit, inv_delays, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([ages_fit; flipud(ages_fit)], [upper_delays; flipud(lower_delays)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds


% Apply moving mean to mylin_curve.y
mylin_curve = myelin_data.myelin; % Adjust window size if needed to smooth

% Define an artificial error range for the fill around mylin_curve
artificial_se = 0.05; % Adjust this value for desired fill width

% Calculate upper and lower bounds for mylin_curve with artificial fill
mylin_curve_upper = myelin_data.upper;
mylin_curve_lower = myelin_data.lower;

% Plot the smoothed line for mylin_curve
plot(myelin_data.Age, mylin_curve, 'b--', 'LineWidth', 1.5);


% Plot the artificial fill for mylin_curve
fill([myelin_data.Age, fliplr(myelin_data.Age)], [mylin_curve_upper, fliplr(mylin_curve_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Age');
ylabel('Myelin');
title('Estimated Myelin (1/\tau^2) vs Myelin Data');
legend('Estimated Myelin', 'Estimator Uncertainty', 'Myelin Data', 'Myelin Uncertainty')
grid on
hold off;

stats.p

%%
clear; clc;

% === Load data ===
dataset = jsondecode(fileread('/mnt/Store/Ronaldo/dev/Data/NewFolder/XIALPHANET.json'));
dataset.Location = '/mnt/Store/Ronaldo/dev/Data/NewFolder';
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

% Z-score outlier removal
z = abs(zscore(delays));
delays = delays(z < 3);
ages = ages(z < 3);

[ages, sortIdx] = sort(ages);
delays = delays(sortIdx);

% === FIGURE 1: Conduction Delay (log fit, normal scale plot) ===
log_delays = log(delays);
X = [ones(size(ages)), ages, ages.^2];
[b1, stats1] = robustfit(X(:,2:end), log_delays);
ages_fit = linspace(min(ages), max(ages), 200)';
X_fit = [ones(size(ages_fit)), ages_fit, ages_fit.^2];
log_y_fit = X_fit * b1;
y_fit = exp(log_y_fit); % back-transform
se_log = sqrt(sum((X_fit * stats1.covb) .* X_fit, 2));
upper = exp(log_y_fit + se_log);
lower = exp(log_y_fit - se_log);

figure('Color','w','Units','normalized','Position',[0.2 0.2 0.6 0.6]);

% Main axis
ax_main = axes('Position', [0.15 0.15 0.65 0.62]); hold on;
fill([ages_fit; flipud(ages_fit)], [upper; flipud(lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ages_fit, y_fit, 'r--', 'LineWidth', 2);

xlabel('Age (Years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conduction Delay (ms)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_main, 'FontSize', 13, 'FontWeight', 'bold'); xlim([5 95]); grid on;

% KDE (top)
ax_top = axes('Position', [0.15 0.79 0.65 0.12]);
[fx, xgrid] = ksdensity(ages);
fill(ax_top, xgrid, fx, [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'k');
axis(ax_top, 'tight'); xlim([5 95]);
set(ax_top, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

% KDE (right)
ax_right = axes('Position', [0.82 0.15 0.12 0.62]);
[fy, ygrid] = ksdensity(delays);
fill(ax_right, fy, ygrid, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'r');
axis(ax_right, 'tight');
set(ax_right, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

% TITLE (annotation ABOVE everything)
annotation('textbox', [0.15, 0.93, 0.7, 0.05], ...
    'String', sprintf('Conduction Delay vs Age (p = %.3g)', stats1.p(3)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');

% === FIGURE 2: Myelin Estimate (log fit, normalized, over myelin data) ===
log_myelin = -2 * log(delays); % log(1/tau^2)
pos = ages< 83;
ages = ages(pos);
log_myelin = log_myelin(pos);

X = [ones(size(ages)), ages, ages.^2];
ages_fit = linspace(0.3, 83, 100)';
X_fit = [ones(size(ages_fit)), ages_fit, ages_fit.^2];
[b2, stats2] = robustfit(X(:,2:end), log_myelin);
log_y_fit = X_fit * b2;
y_fit = exp(log_y_fit); % back-transform
se_log = sqrt(sum((X_fit * stats2.covb) .* X_fit, 2));
upper = exp(log_y_fit + se_log);
lower = exp(log_y_fit - se_log);

% Normalize
y_norm = (y_fit - min(lower)) / (max(upper) - min(lower));
upper = (upper - min(lower)) / (max(upper) - min(lower));
lower = (lower - min(lower)) / (max(upper) - min(lower));

figure('Color','w','Units','normalized','Position',[0.2 0.2 0.6 0.6]);

% Main
ax_main = axes('Position', [0.15 0.15 0.65 0.62]); hold on;
fill([ages_fit; flipud(ages_fit)], [upper; flipud(lower)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ages_fit, y_norm, 'r--', 'LineWidth', 2);

plot(myelin_data.Age, myelin_data.myelin, 'b--', 'LineWidth', 2);
fill([myelin_data.Age, fliplr(myelin_data.Age)], ...
     [myelin_data.upper, fliplr(myelin_data.lower)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Age (Years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Estimated Myelin (1/\tau^2)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_main, 'FontSize', 13, 'FontWeight', 'bold'); xlim([5 95]); grid on;

% KDE (top)
ax_top = axes('Position', [0.15 0.79 0.65 0.12]);
[fx, xgrid] = ksdensity(ages);
fill(ax_top, xgrid, fx, [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'k');
axis(ax_top, 'tight'); xlim([5 95]);
set(ax_top, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

% KDE (right)
ax_right = axes('Position', [0.82 0.15 0.12 0.62]);
[fy, ygrid] = ksdensity(y_norm);
fill(ax_right, fy, ygrid, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'r');
axis(ax_right, 'tight'); ylim([0 1]);
set(ax_right, 'XTick', [], 'YTick', [], 'XColor','none', 'YColor','none'); box off;

% TITLE
annotation('textbox', [0.15, 0.93, 0.7, 0.05], ...
    'String', sprintf('Estimated Myelin vs Myelin Data (p = %.3g)', stats2.p(3)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');

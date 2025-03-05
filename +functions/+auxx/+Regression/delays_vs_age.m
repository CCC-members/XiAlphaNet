% Outlier Elimination and Quadratic Regression Script
% This script reads .mat files from a specified folder, extracts conduction delay and age,
% handles NaN values, eliminates outliers, performs quadratic regression, and plots the results.

% --------------------------- Setup ----------------------------------

% Clear workspace and command window for a fresh start
clear; clc;

%Directory containing .mat files
dataset = jsondecode(fileread('/home/ronaldo/Documents/dev/Data/Results/XIALPHANET.json'));
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
        Mod_Weights = load(fullfile(dataset.Location,participant.SubID,Part_Info.Mod_Weights));
        if participant_age <=15
            delays(index) =  11 * Mod_Weights.Mod_Weights(1);
        else
            delays(index) = 9.5 * Mod_Weights.Mod_Weights(1);
        end
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
        zThreshold = 2.5;

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
        lowerBound = Q1 - 1.5 * IQR;
        upperBound = Q3 + 1.5 * IQR;

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
ages_fit = linspace(0, 100, 100)'; % Same range as LOESS fit
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

% --------------------------- Informative Message -------------------

% Display a message with the number of data points used
fprintf('Quadratic regression performed on %d valid data points after outlier removal.\n', length(cleanDelays));
%%
% Define the ages to fit
% Construct design matrix for quadratic robust fit
invsquare_delays = log(1./(delays.^2));
ages = ages;
X = [ones(size(ages)), ages, ages.^2];
%
 %Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), invsquare_delays); % Exclude intercept from X in robustfit

% Visualize robust regression results

% Evaluate the robust quadratic fit on the same range as ages_fit
ages_fit = linspace(0.3, 82, 100)'; % Same range as LOESS fit
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


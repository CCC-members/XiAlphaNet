%Delayalpha
% --------------------------- Setup ----------------------------------

% Clear workspace and command window for a fresh start
clear; clc;


%Directory containing .mat files
dataset = jsondecode(fileread('D:\data\data\Results\XIALPHANET.json'));
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
        Alpha_estimate = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
        Xi_estimate = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
        %
        threshold_alpha = set_threshold_em(Alpha_estimate.Power);
        pos_alpha = (Alpha_estimate.Power>threshold_alpha);
        alpha_powers(i) =real( mean(Alpha_estimate.Power(pos_alpha)));
        alpha_pfs(i) = real(mean(Alpha_estimate.PAF(pos_alpha)));
        threshold_xi = set_threshold_em(Xi_estimate.Power);
        pos_xi = (Xi_estimate.Power>threshold_xi);
        xi_powers(i) = real(mean(Xi_estimate.Power(pos_xi)));
        if participant_age <=15
            delays(i) =  11 * Mod_Weights.Mod_Weights(1);
        else
            delays(i) = 9.5 * Mod_Weights.Mod_Weights(1);
        end
        index = index +1;
    end
end



%%
delays = delays(:);
ages = ages(:);
alpha_powers = alpha_powers(:);
alpha_pfs = alpha_pfs(:);
xi_powers = xi_powers(:);
%% Alpha_vs_Delays
% Apply outlier removal
clean_ap = alpha_powers;
clean_delays = delays;
% --------------------------- Quadratic Regression ----------------

% Sort data by age for proper plotting
[delays, sortIdx] = sort(clean_delays);
clean_ap = clean_ap(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_ap); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
delays_fit = linspace(min(delays), max(delays), 100)'; % Same range as LOESS fit
robust_fit_interp = (b(1) + b(2) * delays_fit);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(delays_fit)),delays_fit];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  1*se_robust_fit;
lower_bound_robust = robust_fit_interp -  1*se_robust_fit;

% Plot the robust quadratic fit in red
plot(delays_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([delays_fit; flipud(delays_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Delays (ms)');
ylabel('Alpha Power (dB)');
title('Alpha Power vs Delays (\tau)');
legend('Estimated Relationship','Estimator Uncertainty');
grid on;
hold off;

% Apply outlier removal
clean_xip = xi_powers;
clean_delays = delays;
% --------------------------- Quadratic Regression ----------------
clean_xip = clean_xip(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_xip); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
delays_fit = linspace(min(delays), max(delays), 100)'; % Same range as LOESS fit
robust_fit_interp = (b(1) + b(2) * delays_fit);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(delays_fit)),delays_fit];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  1*se_robust_fit;
lower_bound_robust = robust_fit_interp -  1*se_robust_fit;

% Plot the robust quadratic fit in red
plot(delays_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([delays_fit; flipud(delays_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Delays (ms)');
ylabel('Xi Power (dB)');
title('Xi Power vs Delays (\tau)');
legend('Estimated Relationship','Estimator Uncertainty');
grid on;
hold off;

%%
% Apply outlier removal
clean_apf = alpha_pfs;
clean_delays = delays;
% --------------------------- Quadratic Regression ----------------
clean_apf = clean_apf(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_apf); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
delays_fit = linspace(min(delays), max(delays), 100)'; % Same range as LOESS fit
robust_fit_interp = (b(1) + b(2) * delays_fit);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(delays_fit)), delays_fit];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  1*se_robust_fit;
lower_bound_robust = robust_fit_interp -  1*se_robust_fit;

% Plot the robust quadratic fit in red
plot(delays_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([delays_fit; flipud(delays_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Delays (ms)');
ylabel('PAF (Hz)');
title('Peak Alpha Frequency (PAF) vs Delays (\tau)');
legend('Estimated Relationship','Estimator Uncertainty');
grid on;
hold off;
%%
%Alpha_vs_Delays
% Apply outlier removal
clean_ap = alpha_powers;
clean_delays = ages;
% --------------------------- Quadratic Regression ----------------

% Sort data by age for proper plotting
[delays, sortIdx] = sort(clean_delays);
clean_ap = clean_ap(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays, delays.^2];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_ap); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
delays_fit = linspace(min(delays), max(delays), 100)'; % Same range as LOESS fit
robust_fit_interp = (b(1) + b(2) * delays_fit+b(3) * delays_fit.^2);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(delays_fit)), delays_fit,delays_fit.^2];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  1*se_robust_fit;
lower_bound_robust = robust_fit_interp -  1*se_robust_fit;

% Plot the robust quadratic fit in red
plot(delays_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([delays_fit; flipud(delays_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Age (Years)');
ylabel('Alpha Power (dB)');
title('Alpha Power vs Age');
legend('Estimated Relationship','Estimator Uncertainty');
grid on;
hold off;

% Apply outlier removal
clean_xip = xi_powers;
clean_delays = ages;
% --------------------------- Quadratic Regression ----------------
clean_xip = clean_xip(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays, delays.^2];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_xip); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
delays_fit = linspace(min(delays), max(delays), 100)'; % Same range as LOESS fit
robust_fit_interp = (b(1) + b(2) * delays_fit+b(3) * delays_fit.^2);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(delays_fit)), delays_fit,delays_fit.^2];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  1*se_robust_fit;
lower_bound_robust = robust_fit_interp -  1*se_robust_fit;

% Plot the robust quadratic fit in red
plot(delays_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([delays_fit; flipud(delays_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Age (Years)');
ylabel('Xi Power (dB)');
title('Xi Power vs Age');
legend('Estimated Relationship','Estimator Uncertainty');
grid on;
hold off;


% Apply outlier removal
clean_apf = alpha_pfs;
clean_delays = ages;
% --------------------------- Quadratic Regression ----------------
clean_apf = clean_apf(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays, delays.^2];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_apf); % Exclude intercept from X in robustfit

% Visualize robust regression results
figure;
%scatter(ages, transformed_delays, 'filled', 'SizeData', 1.5);
hold on;

% Evaluate the robust quadratic fit on the same range as ages_fit
delays_fit = linspace(min(delays), max(delays), 100)'; % Same range as LOESS fit
robust_fit_interp = (b(1) + b(2) * delays_fit+b(3) * delays_fit.^2);

% Calculate standard error for each fitted value
cov_b = stats.covb; % Covariance matrix of coefficients
X_fit = [ones(size(delays_fit)), delays_fit,delays_fit.^2];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  1*se_robust_fit;
lower_bound_robust = robust_fit_interp -  1*se_robust_fit;

% Plot the robust quadratic fit in red
plot(delays_fit, robust_fit_interp, 'r--', 'LineWidth', 1.5); % Robust quadratic fit line

% Fill the area between the error bounds for the robust fit with red color
fill([delays_fit; flipud(delays_fit)], [upper_bound_robust; flipud(lower_bound_robust)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded red area between bounds
%scatter(ages,delays,'blue','filled')
hold off
% --------------------------- Plotting -------------------------------

% Enhance plot aesthetics
xlabel('Age (Years)');
ylabel('PAF (Hz)');
title('Peak Alpha Frequency (PAF) vs Age');
legend('Estimated Relationship','Estimator Uncertainty');
grid on;
hold off;
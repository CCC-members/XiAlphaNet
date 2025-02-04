%Delayalpha
% --------------------------- Setup ----------------------------------

% Clear workspace and command window for a fresh start
clear; clc;

% Directory containing .mat files
dataFolder = 'Data/Model_Parameters/Control2/'; % Replace with your actual folder path

% Initialize arrays for storing delays and ages
delays = [];
ages = [];
alpha_powers= [];
xi_powers= [];
alpha_pfs = [];
% Get a list of all .mat files in the folder
files = dir(fullfile(dataFolder, '*.mat'));

% Check if any .mat files are found
if isempty(files)
    error('No .mat files found in the specified folder: %s', dataFolder);
end

% --------------------------- Data Extraction -----------------------

% Loop through each file and load the data
perm1 = randperm(length(files),length(files))
for i = perm1
    % Load the .mat file
    filePath = fullfile(dataFolder, files(i).name);
    data = load(filePath);
    
    % Extract conduction delay and age
    if isfield(data, 'x') && isfield(data.x, 'Lambda_DC') && isfield(data.x, 'Age')
        delay = 9.6*data.x.Lambda_DC(1);%data.d.Lambda_DC(1);
        age = data.x.Age;
        [e,a,s2] = x2v(data.x.Solution);
        pos_alpha = (a(:,1)>0); 
        alpha_power =real( mean(a(pos_alpha,1)));
        alpha_pf = real(mean(a(pos_alpha,4)));
        pos_xi = (e(:,1)>0); 
        xi_power = real(mean(e(pos_alpha,1)));

        % Check for NaN values and ensure they are numeric
        if isnumeric(delay) && isnumeric(age) && ~isnan(delay) && ~isnan(age)
            delays(end+1) = delay; %#ok<SAGROW>
            ages(end+1) = age; %#ok<SAGROW>
            alpha_powers(end+1) = alpha_power;
            alpha_pfs(end+1) = alpha_pf;
            xi_powers(end+1) = xi_power;
        else
            fprintf('NaN or non-numeric value found in file: %s. Skipping this entry.\n', files(i).name);
        end
    else
        fprintf('Missing fields in file: %s. Skipping this file.\n', files(i).name);
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
X_fit = [ones(size(delays_fit)), delays_fit];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  0.5*se_robust_fit;
lower_bound_robust = robust_fit_interp -  0.5*se_robust_fit;

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
%%
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
X_fit = [ones(size(delays_fit)), delays_fit];
se_robust_fit = sqrt(sum((X_fit * cov_b) .* X_fit, 2));

% Calculate upper and lower bounds for the robust fit curve
upper_bound_robust = robust_fit_interp +  0.5*se_robust_fit;
lower_bound_robust = robust_fit_interp -  0.5*se_robust_fit;

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
upper_bound_robust = robust_fit_interp +  0.5*se_robust_fit;
lower_bound_robust = robust_fit_interp -  0.5*se_robust_fit;

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
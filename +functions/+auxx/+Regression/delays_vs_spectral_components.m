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
delays0 = delays(:);
ages0 = ages(:);
alpha_powers0 = alpha_powers(:);
alpha_pfs0 = alpha_pfs(:);
xi_powers0 = xi_powers(:);

delays = delays0(:);
ages = ages0(:);
alpha_powers = alpha_powers0(:);
alpha_pfs = alpha_pfs0(:);
xi_powers = xi_powers0(:);
%% Alpha_vs_Delays
% Apply outlier removal
clean_ap = (alpha_powers-min(alpha_powers))./(max(alpha_powers)-min(alpha_powers));
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
clean_apf = (alpha_pfs-min(alpha_pfs))./(max(alpha_pfs)-min(alpha_pfs));
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SINGLE SCRIPT: ROBUST STATISTICAL ANALYSIS WITH MIN–MAX SCALING       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Data Setup (Assumes the following exist in workspace):
%   delays       = column vector (Nx1)
%   alpha_powers = column vector (Nx1)
%   xi_powers    = column vector (Nx1)
%   alpha_pfs    = column vector (Nx1)
%
% Ensure each is a column vector:
delays       = delays(:);
alpha_powers = alpha_powers(:);
xi_powers    = xi_powers(:);
alpha_pfs    = alpha_pfs(:);

% Set random seed for reproducible bootstrap/permutation results (optional)
rng(1234);

% Define some parameters for bootstrapping & permutation tests:
numBoots   = 1000;   % number of bootstrap samples
alphaLevel = 0.05;   % for 95% CI
numPerms   = 10000;  % for permutation test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             A) ALPHA POWER vs. DELAYS  (Theil–Sen + Bootstraps)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('=== ALPHA POWER vs. DELAYS ===');

% 1) Outlier Removal (IQR-based) on alpha_powers (y)
x = delays;
y = alpha_powers;

% Remove NaNs just in case
validIdx = ~isnan(x) & ~isnan(y);
x = x(validIdx);
y = y(validIdx);

% IQR-based outliers on y
Q1 = quantile(y, 0.25);
Q3 = quantile(y, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5 * IQR_val;
upperBound = Q3 + 1.5 * IQR_val;

idxOutliers = (y < lowerBound) | (y > upperBound);
xClean = x(~idxOutliers);
yClean = y(~idxOutliers);

fprintf('Removed %d outliers (alpha_powers) via IQR.\n', sum(idxOutliers));

% 2) Min–max scale xClean and yClean to [0,1]
xMin = min(xClean);  xMax = max(xClean);
yMin = min(yClean);  yMax = max(yClean);

% Avoid dividing by zero if data are constant
if xMax > xMin
    xClean = (xClean - xMin) / (xMax - xMin);
else
    xClean = zeros(size(xClean));
end

if yMax > yMin
    yClean = (yClean - yMin) / (yMax - yMin);
else
    yClean = zeros(size(yClean));
end

% 3) Theil–Sen slope & intercept on scaled data
n = length(xClean);
allSlopes = [];
idxSlope = 0;
for i = 1:n-1
    for j = i+1:n
        if xClean(j) ~= xClean(i)
            idxSlope = idxSlope + 1;
            allSlopes(idxSlope) = (yClean(j) - yClean(i)) / ...
                                  (xClean(j) - xClean(i));
        end
    end
end
theilSenSlope = median(allSlopes);

% Intercept = median of (y - slope*x)
allIntercepts = yClean - theilSenSlope*xClean;
theilSenIntercept = median(allIntercepts);

fprintf('Theil-Sen Slope = %.6f\n', theilSenSlope);
fprintf('Theil-Sen Intercept = %.6f\n', theilSenIntercept);

% 4) Bootstrap CI for Theil–Sen slope
bootSlopes = zeros(numBoots,1);
for b = 1:numBoots
    idxB = randi(n, [n,1]);  % sample with replacement
    xB = xClean(idxB);
    yB = yClean(idxB);

    localSlopes = [];
    idxBSlope = 0;
    for i2 = 1:n-1
        for j2 = i2+1:n
            if xB(j2) ~= xB(i2)
                idxBSlope = idxBSlope + 1;
                localSlopes(idxBSlope) = (yB(j2) - yB(i2)) / ...
                                         (xB(j2) - xB(i2));
            end
        end
    end
    bootSlopes(b) = median(localSlopes);
end

ciLow = prctile(bootSlopes, 100*(alphaLevel/2));
ciHigh= prctile(bootSlopes, 100*(1 - alphaLevel/2));
fprintf('95%% Bootstrap CI for Theil-Sen slope: [%.6f, %.6f]\n', ciLow, ciHigh);

if (ciLow>0 && ciHigh>0) || (ciLow<0 && ciHigh<0)
    disp('Slope is significantly different from zero (does not cross 0).');
else
    disp('Slope is NOT significantly different from zero (CI crosses 0).');
end

% 5) Spearman Correlation on scaled xClean and yClean
[spearmanRho, spearmanP] = corr(xClean, yClean, 'Type','Spearman');
fprintf('Spearman rho = %.4f, p = %.4g\n', spearmanRho, spearmanP);

% 6) Permutation Test for correlation (two-sided)
observedRho = spearmanRho;
permVals = zeros(numPerms,1);
for pIdx = 1:numPerms
    yPerm = yClean(randperm(n));
    permVals(pIdx) = corr(xClean, yPerm, 'Type','Spearman');
end
pPerm = mean(abs(permVals) >= abs(observedRho));
fprintf('Permutation test (two-sided) for Spearman corr: p = %.4f\n\n', pPerm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             B) XI POWER vs. DELAYS  (Same Steps, Scaled)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('=== XI POWER vs. DELAYS ===');

x = delays;
y = xi_powers;

validIdx = ~isnan(x) & ~isnan(y);
x = x(validIdx);
y = y(validIdx);

Q1 = quantile(y, 0.25);
Q3 = quantile(y, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5 * IQR_val;
upperBound = Q3 + 1.5 * IQR_val;

idxOutliers = (y < lowerBound) | (y > upperBound);
xClean = x(~idxOutliers);
yClean = y(~idxOutliers);

fprintf('Removed %d outliers (xi_powers) via IQR.\n', sum(idxOutliers));

% Scale xClean, yClean
xMin = min(xClean);  xMax = max(xClean);
yMin = min(yClean);  yMax = max(yClean);

if xMax > xMin
    xClean = (xClean - xMin) / (xMax - xMin);
else
    xClean = zeros(size(xClean));
end

if yMax > yMin
    yClean = (yClean - yMin) / (yMax - yMin);
else
    yClean = zeros(size(yClean));
end

n = length(xClean);
allSlopes = [];
idxSlope = 0;
for i = 1:n-1
    for j = i+1:n
        if xClean(j) ~= xClean(i)
            idxSlope = idxSlope + 1;
            allSlopes(idxSlope) = (yClean(j) - yClean(i)) / ...
                                  (xClean(j) - xClean(i));
        end
    end
end
theilSenSlope = median(allSlopes);
allIntercepts = yClean - theilSenSlope*xClean;
theilSenIntercept = median(allIntercepts);

fprintf('Theil-Sen Slope = %.6f\n', theilSenSlope);
fprintf('Theil-Sen Intercept = %.6f\n', theilSenIntercept);

bootSlopes = zeros(numBoots,1);
for b = 1:numBoots
    idxB = randi(n, [n,1]);
    xB = xClean(idxB);
    yB = yClean(idxB);

    localSlopes = [];
    idxBSlope = 0;
    for i2 = 1:n-1
        for j2 = i2+1:n
            if xB(j2) ~= xB(i2)
                idxBSlope = idxBSlope + 1;
                localSlopes(idxBSlope) = (yB(j2) - yB(i2)) / ...
                                         (xB(j2) - xB(i2));
            end
        end
    end
    bootSlopes(b) = median(localSlopes);
end

ciLow = prctile(bootSlopes, 100*(alphaLevel/2));
ciHigh= prctile(bootSlopes, 100*(1 - alphaLevel/2));
fprintf('95%% Bootstrap CI for Theil-Sen slope: [%.6f, %.6f]\n', ciLow, ciHigh);

if (ciLow>0 && ciHigh>0) || (ciLow<0 && ciHigh<0)
    disp('Slope is significantly different from zero.');
else
    disp('Slope is NOT significantly different from zero.');
end

[spearmanRho, spearmanP] = corr(xClean, yClean, 'Type','Spearman');
fprintf('Spearman rho = %.4f, p = %.4g\n', spearmanRho, spearmanP);

observedRho = spearmanRho;
permVals = zeros(numPerms,1);
for pIdx = 1:numPerms
    yPerm = yClean(randperm(n));
    permVals(pIdx) = corr(xClean, yPerm, 'Type','Spearman');
end
pPerm = mean(abs(permVals) >= abs(observedRho));
fprintf('Permutation test (two-sided): p = %.4f\n\n', pPerm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            C) PEAK ALPHA FREQUENCY vs. DELAYS  (Scaled)                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('=== PAF vs. DELAYS ===');

x = delays;
y = alpha_pfs;

validIdx = ~isnan(x) & ~isnan(y);
x = x(validIdx);
y = y(validIdx);

Q1 = quantile(y, 0.25);
Q3 = quantile(y, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5 * IQR_val;
upperBound = Q3 + 1.5 * IQR_val;

idxOutliers = (y < lowerBound) | (y > upperBound);
xClean = x(~idxOutliers);
yClean = y(~idxOutliers);

fprintf('Removed %d outliers (paf) via IQR.\n', sum(idxOutliers));

% Scale xClean, yClean
xMin = min(xClean);  xMax = max(xClean);
yMin = min(yClean);  yMax = max(yClean);

if xMax > xMin
    xClean = (xClean - xMin) / (xMax - xMin);
else
    xClean = zeros(size(xClean));
end

if yMax > yMin
    yClean = (yClean - yMin) / (yMax - yMin);
else
    yClean = zeros(size(yClean));
end

n = length(xClean);
allSlopes = [];
idxSlope = 0;
for i = 1:n-1
    for j = i+1:n
        if xClean(j) ~= xClean(i)
            idxSlope = idxSlope + 1;
            allSlopes(idxSlope) = (yClean(j) - yClean(i)) / ...
                                  (xClean(j) - xClean(i));
        end
    end
end
theilSenSlope = median(allSlopes);
allIntercepts = yClean - theilSenSlope*xClean;
theilSenIntercept = median(allIntercepts);

fprintf('Theil-Sen Slope = %.6f\n', theilSenSlope);
fprintf('Theil-Sen Intercept = %.6f\n', theilSenIntercept);

bootSlopes = zeros(numBoots,1);
for b = 1:numBoots
    idxB = randi(n, [n,1]);
    xB = xClean(idxB);
    yB = yClean(idxB);

    localSlopes = [];
    idxBSlope = 0;
    for i2 = 1:n-1
        for j2 = i2+1:n
            if xB(j2) ~= xB(i2)
                idxBSlope = idxBSlope + 1;
                localSlopes(idxBSlope) = (yB(j2) - yB(i2)) / ...
                                         (xB(j2) - xB(i2));
            end
        end
    end
    bootSlopes(b) = median(localSlopes);
end

ciLow = prctile(bootSlopes, 100*(alphaLevel/2));
ciHigh= prctile(bootSlopes, 100*(1 - alphaLevel/2));
fprintf('95%% Bootstrap CI for Theil-Sen slope: [%.6f, %.6f]\n', ciLow, ciHigh);

if (ciLow>0 && ciHigh>0) || (ciLow<0 && ciHigh<0)
    disp('Slope is significantly different from zero.');
else
    disp('Slope is NOT significantly different from zero.');
end

[spearmanRho, spearmanP] = corr(xClean, yClean, 'Type','Spearman');
fprintf('Spearman rho = %.4f, p = %.4g\n', spearmanRho, spearmanP);

observedRho = spearmanRho;
permVals = zeros(numPerms,1);
for pIdx = 1:numPerms
    yPerm = yClean(randperm(n));
    permVals(pIdx) = corr(xClean, yPerm, 'Type','Spearman');
end
pPerm = mean(abs(permVals) >= abs(observedRho));
fprintf('Permutation test (two-sided): p = %.4f\n', pPerm);

disp('=== END OF SCRIPT ===');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MULTIPLE LINEAR (ROBUST) MODEL: PAF = β0 + β1*(Delays) + β2*(Age)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Suppose you already have three column vectors in your workspace:
%   paf    = Nx1 vector (Peak Alpha Frequency measurements)
%   delays = Nx1 vector (Delays in ms or some unit)
%   age    = Nx1 vector (Age in years)

% 1) Ensure they are column vectors
paf    =(alpha_pfs(:)-min(alpha_pfs(:)))./(max(alpha_pfs(:))-min(alpha_pfs(:)));
delays = 1./delays0(:).^1;
delays = (delays-min(delays))./(max(delays)-min(delays));
age    = ages(:);

% 2) (Optional) Remove outliers or NaNs if desired:
validIdx = ~isnan(paf) & ~isnan(delays) & ~isnan(age);
pafC    = paf(validIdx);
delaysC = delays(validIdx);
ageC    = age(validIdx);

% Example IQR-based outlier removal on PAF:
Q1 = quantile(pafC, 0.25);
Q3 = quantile(pafC, 0.75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - 1.5*IQR_val;
upperBound = Q3 + 1.5*IQR_val;

idxOut = (pafC < lowerBound) | (pafC > upperBound);
fprintf('Removed %d outliers (PAF) via IQR.\n', sum(idxOut));

pafClean    = pafC(~idxOut);
delaysClean = delaysC(~idxOut);
ageClean    = ageC(~idxOut);


% ===================== 1) CONSTRUCT NONLINEAR DESIGN MATRIX ==============
% We want columns for Age, Age^2, Delays, Delays^2, and Age*Delays.
% robustfit automatically adds an intercept, so we do NOT include a column of ones.

X = [ ageClean, (ageClean.^2), ...
      delaysClean,delaysClean.^2,ageClean.*delaysClean];

% ===================== 2) FIT ROBUST MODEL ===============================
[b, stats] = robustfit(X, pafClean);

% b(1) = intercept
% b(2) = coefficient for Age
% b(3) = coefficient for Age^2
% b(4) = coefficient for Delays
% b(5) = coefficient for Delays^2
% b(6) = coefficient for Age*Delays

% stats structure includes p-values, residuals, standard errors, etc.
pVals = stats.p;    % p-values
SE    = stats.se;   % standard errors

% ===================== 3) DISPLAY RESULTS ================================
fprintf('\n=== Nonlinear Robust Model: PAF ~ Age + Age^2 + Delays ===\n');
fprintf('Intercept        = %.4f (SE=%.4f, p=%.4g)\n', b(1), SE(1), pVals(1));
fprintf('Age              = %.4f (SE=%.4f, p=%.4g)\n', b(2), SE(2), pVals(2));
fprintf('Age^2            = %.4f (SE=%.4f, p=%.4g)\n', b(3), SE(3), pVals(3));
fprintf('Delays           = %.4f (SE=%.4f, p=%.4g)\n', b(4), SE(4), pVals(4));
fprintf('Delays^2         = %.4f (SE=%.4f, p=%.4g)\n', b(5), SE(5), pVals(5));

% Interpretation:
%  - If pVals(2) < 0.05 => Age is a significant linear term
%  - pVals(3) < 0.05 => Age^2 is significant (curvature effect)
%  - pVals(4) < 0.05 => Delays is significant
%  - pVals(5) < 0.05 => Delays^2 is significant (curvature effect)
%  - pVals(6) < 0.05 => The interaction of Age & Delays is significant

% ===================== 4) OPTIONAL: DIAGNOSTICS / PLOTS ==================
% e.g., check residual distribution, or generate a 3D surface plot if you want.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NONLINEAR ROBUST MODEL + QUADRATIC PLOTS FOR PAF VS AGE & DELAYS VS AGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Steps 0 through 4 from your script remain the same, up to "disp('Done.');"

% (We assume you already have variables: ageClean, pafClean, delaysClean)
% which are scaled between 0-1 and outliers removed, as in your script.

%% 5) SCATTER PLOT OF THE DATA
figure; hold on;
plot(ageClean, pafClean, '.', 'Color', 'r');
plot(ageClean, delaysClean, '.', 'Color', 'b');

% 6) QUADRATIC FIT FOR PAF vs. AGE
% We want a polynomial pPAF of degree 2 such that:
%   pafClean ~= pPAF(1)*age^2 + pPAF(2)*age + pPAF(3)
pPAF = polyfit(ageClean, pafClean, 2);

% Create a smooth age grid for plotting the fitted curve
ageFit = linspace(min(ageClean), max(ageClean), 200);

% Evaluate the polynomial on this grid
pafFit = polyval(pPAF, ageFit);

% Plot the PAF-vs-Age quadratic in red
plot(ageFit, pafFit, 'r-', 'LineWidth', 4);

% 7) QUADRATIC FIT FOR DELAYS vs. AGE
% Similarly for delays
pDel = polyfit(ageClean, delaysClean, 2);
delFit = polyval(pDel, ageFit);

% Plot the Delays-vs-Age quadratic in blue
plot(ageFit, delFit, 'b-', 'LineWidth', 4);

% 8) FINALIZE THE PLOT
xlabel('Age');
ylabel('Scaled PAF (red) or Scaled Delays (blue)');
title('Quadratic Fits for PAF and Delays vs. Age');
legend({'PAF data','Delays data','PAF Quadratic Fit','Delays Quadratic Fit'}, ...
       'Location','best');
grid on;
hold off;


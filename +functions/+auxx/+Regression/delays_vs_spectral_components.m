%Delayalpha
% --------------------------- Setup ----------------------------------

% Clear workspace and command window for a fresh start
clear; clc;


%Directory containing .mat files
dataset = jsondecode(fileread('/Users/ronald/Downloads/Results/XIALPHANET.json'));
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
for i=1:length(dataset.Participants)
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
        threshold_alpha =  prctile(Alpha_estimate.Power,90);%set_threshold_em(Alpha_estimate.Power);
        pos_alpha = (Alpha_estimate.Power>threshold_alpha);
        alpha_powers(i) =real( mean(Alpha_estimate.Power(pos_alpha)));
        alpha_pfs(i) = real(mean(Alpha_estimate.PAF(pos_alpha)));
        threshold_xi = prctile(Xi_estimate.Power,90);%set_threshold_em(Xi_estimate.Power);
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
%%
delays = delays0(:);
ages = ages0(:);
alpha_powers = alpha_powers0(:);
alpha_pfs = alpha_pfs0(:);
xi_powers = xi_powers0(:);
%% Alpha_vs_Delays
% Apply outlier removal
clean_ap = log(10^(-6)+alpha_powers);
clean_ap = (clean_ap-min(clean_ap))./(max(clean_ap)-min(clean_ap));
clean_delays = log(delays0.^(-2));
clean_delays = (clean_delays-min(clean_delays))./(max(clean_delays)-min(clean_delays));
% --------------------------- Quadratic Regression ----------------

% Sort data by age for proper plotting
[delays, sortIdx] = sort(clean_delays);
clean_ap = clean_ap(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_ap); % Exclude intercept from X in robustfit
stats.p
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
clean_xip = log(10^(-6)+xi_powers);
clean_xip = (clean_xip-min(clean_xip))./(max(clean_xip)-min(clean_xip));
clean_delays = log(delays0.^(-2));
clean_delays = (clean_delays-min(clean_delays))./(max(clean_delays)-min(clean_delays));
% --------------------------- Quadratic Regression ----------------
clean_xip = clean_xip(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_xip); % Exclude intercept from X in robustfit
stats.p
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

% Apply outlier removal
clean_apf = log(10^(-6)+alpha_pfs);
clean_apf = (clean_apf-min(clean_apf))./(max(clean_apf)-min(clean_apf));
clean_delays = log(delays0.^(-2));
clean_delays = (clean_delays-min(clean_delays))./(max(clean_delays)-min(clean_delays));
% --------------------------- Quadratic Regression ----------------

% Sort data by age for proper plotting
[delays, sortIdx] = sort(clean_delays);
clean_apf = clean_apf(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_apf); % Exclude intercept from X in robustfit
stats.p
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
%Alpha_vs_Age
% Apply outlier removal
clean_ap = log(1+alpha_powers);
clean_delays = ages;
% --------------------------- Quadratic Regression ----------------

% Sort data by age for proper plotting
[delays, sortIdx] = sort(clean_delays);
clean_ap = clean_ap(sortIdx);

% Construct design matrix for quadratic robust fit
X = [ones(size(delays)), delays, delays.^2];

% Perform robust regression on the transformed data with quadratic model
[b, stats] = robustfit(X(:, 2:end), clean_ap); % Exclude intercept from X in robustfit
stats.p
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
delays = 1./delays0(:).^2;
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
ageFit = linspace(0, 80, 200);

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
%%
% 0) Ensure columns, remove NaNs if needed
paf    = alpha_pfs(:);
delays = delays0(:);
age    = ages0(:);

validIdx = ~isnan(paf) & ~isnan(delays) & ~isnan(age);
pafClean    = paf(validIdx);
delaysClean = delays(validIdx);
ageClean    = age(validIdx);

% 1) Compute 1/delay^2
invDelays2 = delaysClean.^1;

fprintf('\n=== PART A: PARTIAL CORRELATION (CONTROLLING FOR AGE) ===\n');

% A1) Partial correlation controlling only for Age (linear)
% partialcorr(X, Y, Z) returns the correlation between X and Y
% after removing the effect of Z.
% We'll use Spearman's rank correlation to be more robust to outliers.
[rPartial_linear, pPartial_linear] = partialcorr(pafClean, invDelays2, ageClean, ...
                                                'type', 'Spearman');
fprintf('Partial Corr (PAF, 1/Del^2 | Age): r=%.4f, p=%.4g\n', ...
        rPartial_linear, pPartial_linear);

% A2) Partial correlation controlling for Age and Age^2
% We'll build a matrix for the control variables [Age, Age^2].
Z = [ageClean, ageClean.^2];
[rPartial_quad, pPartial_quad] = partialcorr(pafClean, invDelays2, Z, ...
                                             'type', 'Spearman');
fprintf('Partial Corr (PAF, 1/Del^2 | Age, Age^2): r=%.4f, p=%.4g\n', ...
        rPartial_quad, pPartial_quad);

% Interpretation:
% If these partial correlations remain significant (p<0.05), it suggests
% that 1/Del^2 correlates with PAF even after factoring out age.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART B: MULTIPLE (ROBUST) REGRESSION: PAF ~ Age + Age^2 + (1/Del^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== PART B: MULTIPLE ROBUST REGRESSION ===\n');

% B1) Construct design matrix for robustfit.
% robustfit automatically includes an intercept, so just columns for:
%   1) Age
%   2) Age^2
%   3) 1/Del^2
X = [ ageClean, ageClean.^2, invDelays2 ];

% B2) Fit robust regression: pafClean = b0 + b1*Age + b2*Age^2 + b3*(1/Del^2)
[b, stats] = robustfit(X, pafClean);

% b(1) = intercept
% b(2) = coefficient for Age
% b(3) = coefficient for Age^2
% b(4) = coefficient for (1/Del^2)

pVals = stats.p;
SEs   = stats.se;
coeffNames = {'(Intercept)','Age','Age^2','1/delay^2'};

fprintf('Coefficients (Robust Fit):\n');
for i = 1:numel(coeffNames)
    fprintf('  %s = %.4f (SE=%.4f, p=%.4g)\n', ...
            coeffNames{i}, b(i), SEs(i), pVals(i));
end

% Interpretation:
%  - If pVals(4) < 0.05 => (1/Del^2) is a significant predictor of PAF
%    even after controlling for Age and Age^2.
%  - Age or Age^2 being significant indicates an inverted-U or other
%    shape in the PAF vs. Age relationship.

disp('=== Analysis Complete ===');

%%
%% --------------------------- Alpha Power vs Delays^-2 ---------------------------
% 1) Log-transform and sort
AP = 10*log10(1e-6 + alpha_powers);          % log(Alpha Power)
%AP =  (AP(:)-min(AP(:)))./(max(AP(:))-min(AP(:)));
Del = (delays0);               % log(tau^-2)
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
title('Alpha Power ~ \tau  (Robust Linear & Quadratic)');
xlim([4 18])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;
% --------------------------- Peak Alpha Frequency vs Delays^-2 ---------------------------
% 1) Log-transform and sort
AP = alpha_pfs; % log(Alpha Power)
%AP =  (AP(:)-min(AP(:)))./(max(AP(:))-min(AP(:)));
Del = (delays0);               % log(tau^-2)
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
title('Peak Alpha Frequency (PAF) ~ \tau  (Robust Linear & Quadratic)');
xlim([4 18])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

% --------------------------- Xi Power vs Delays^-2 ---------------------------
XIP = 10*log10(10^(-6)+xi_powers); 
%XIP =  (XIP(:)-min(XIP(:)))./(max(XIP(:))-min(XIP(:)));
Del = (delays0);               % log(tau^-2)
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
title('Xi Power ~ \tau (Robust Linear & Quadratic)');
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
title('Alpha Power ~ age  (Robust Linear & Quadratic)');
xlim([5 95])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;
% --------------------------- Peak Alpha Frequency vs Delays^-2 ---------------------------
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
title('Peak Alpha Frequency (PAF) ~ age  (Robust Linear & Quadratic)');
xlim([5 95])
% Legend with parameter values, RMSE, BIC and p-values
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

% --------------------------- Xi Power vs Delays^-2 ---------------------------
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
title('Xi Power ~ age (Robust Linear & Quadratic)');
xlim([5 95])
legLin = sprintf('Linear: Intercept=%.3f (p=%.3g), Slope=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_lin(1), stats_lin.p(1), b_lin(2), stats_lin.p(2), rmse_lin, bic_lin);
legQuad = sprintf('Quadratic: b0=%.3f (p=%.3g), b1=%.3f (p=%.3g), b2=%.3f (p=%.3g), RMSE=%.3f, BIC=%.3f', ...
    b_quad(1), stats_quad.p(1), b_quad(2), stats_quad.p(2), b_quad(3), stats_quad.p(3), rmse_quad, bic_quad);

legend({legLin, 'Linear \pm1SE', legQuad, 'Quadratic \pm1SE'}, ...
    'Location', 'best');
grid on; hold off;

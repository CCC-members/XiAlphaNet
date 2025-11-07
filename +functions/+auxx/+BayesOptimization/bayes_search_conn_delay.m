function lambda_opt_dc = bayes_search_conn_delay(lambda_space, Ne, Nr, Nw, freq, Cross, ...
    K, D, C, indx_parallel, BayesIter_Delay, x0, Lipschitz)
% BAYES_SEARCH_CONN_DELAY
% Bayesian search for (lambda_D, lambda_C) scaling parameters.
%
% - Deterministic objective (evaluateF)
% - Log-normal priors on lambda1 (delay scaling) and lambda2 (connectivity scaling)
% - Optional parallel execution
%
% INPUTS:
%   lambda_space   : [2 x 2] bounds for [lambda1; lambda2]
%   Ne, Nr, Nw     : model dimensions
%   freq           : frequency grid
%   Cross          : cross-spectral data
%   K, D, C        : model operators (kernel, delays, connectivity)
%   indx_parallel  : flag for parallel BO
%   BayesIter_Delay: number of BO evaluations
%   x0, Lipschitz  : parameters for evaluateF
%
% OUTPUT:
%   lambda_opt_dc  : [lambda1, lambda2] optimal values

% Use all frequencies
Nsfreq = length(freq);

% Objective function
objectiveFunc = @(lambda) modelObjective(lambda, ...
    Ne, Nr, Nw, freq, Cross, K, D, C, x0, Lipschitz, Nsfreq, lambda_space);

% Define domains (log-scale search)
ld_domain = optimizableVariable('l1', ...
    [lambda_space(1,1), lambda_space(1,2)], 'Transform', 'log');
lc_domain = optimizableVariable('l2', ...
    [lambda_space(2,1), lambda_space(2,2)], 'Transform', 'log');

% Configure Bayesian optimization
results = bayesopt(objectiveFunc, [ld_domain, lc_domain], ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'MaxObjectiveEvaluations', BayesIter_Delay, ...
    'IsObjectiveDeterministic', true, ...
    'NumSeedPoints', 10, ...
    'UseParallel', indx_parallel, ...
    'Verbose', 0, ...
    'PlotFcn', []);

% Extract optimum
lambda_opt_dc = [results.XAtMinObjective.l1, results.XAtMinObjective.l2];
end


% ======================= Objective Function =======================
function Fval = modelObjective(lambda, Ne, Nr, Nw, freq, Cross, ...
    K, D, C, x0, Lipschitz, Nsfreq, lambda_space)

% Keep full imports so dependencies are visible to workers
import functions.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.Regularization.*
import functions.auxx.TOperator.*
import functions.StochasticFISTA.*
import functions.auxx.RegSpace.*
import functions.auxx.StochasticEval.*
import functions.FunctionGrandientProx.*

% Frequency subsampling
[nsf_band, sw, sp] = sample_frequencies(freq, 0, Nsfreq);

% Extract scaling parameters
lambda1 = lambda.l1;  % delay scaling
lambda2 = lambda.l2;  % connectivity scaling

% Build compact model
params.Dimensions.Ne   = Ne;
params.Dimensions.Nr   = Nr;
params.Dimensions.Nv   = Nr;
params.Dimensions.Nw   = Nw;
params.Data.Cross      = Cross;
params.Data.freq       = freq;
params.Compact_Model.K = K;
params.Compact_Model.D = lambda1 * D;
params.Compact_Model.C = lambda2 * C;
params.Parallel.T      = 0;

% Evaluate operator
T = Teval(params);

% Deterministic evaluation
FS = evaluateF(x0, Ne, T, sw, sp, nsf_band, Cross);

% Typical scale for penalties
F_scale = abs(FS);
beta    = 0.01;  % relative penalty strength

% Penalty scaling for lambda1
delta_max1 = max((log(lambda_space(1,1)))^2, (log(lambda_space(1,2)))^2);
alpha1     = (beta * F_scale) / delta_max1;

% Penalty scaling for lambda2
eps_val    = 1e-6;
delta_max2 = max((log(lambda_space(2,1) + eps_val))^2, ...
    (log(lambda_space(2,2)))^2);
alpha2     = (beta * F_scale) / delta_max2;

% Log-normal priors (centered at 1)
log_diff1 = log(lambda1);
log_diff2 = log(lambda2);
penalty   = alpha1 * (log_diff1.^2) + alpha2 * (log_diff2.^2);

% Augmented objective
Fval = FS + penalty;
end

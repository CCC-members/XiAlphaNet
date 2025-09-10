function [lambda_opt_dc] = bayes_search_conn_delay(lambda_space, Ne, Nr, Nw, freq, Cross, BayesIter_Conn, K, D, C, indx_parallel, BayesIter_Delay, x0, Lipschitz)
% Variance-reduced Bayesian search for (lambda_D, lambda_C)
% - Common Random Numbers (fixed RNG seeds) inside the objective
% - R replicates with trimmed-mean aggregation (robust to outliers)
% - Deterministic objective to the BO GP ('IsObjectiveDeterministic' = true)
% - Optional parallelism (safe, since objective is deterministic)
% - Final deterministic recheck of BO's top-k candidates with full data

% ---------- Minimal, tunable controls (no new args added) ----------
R_reps        = 50;        % # replicates per lambda for robust averaging
trimFrac      = 0.15;     % 15% trimmed mean (total trim; split on both tails)
baseSeed      = 123456;   % fixed base seed for CRN
NsfreqFrac    = 0.9;     % fraction of frequencies used during BO
TopK_recheck  = 7;        % reevaluate this many best-BO candidates at the end
R_recheck     = 9;        % replicates during the deterministic recheck
% -------------------------------------------------------------------

% Precompute CRN replicate seeds (shared across ALL lambda evals)
repSeeds = baseSeed + (1:R_reps);

% Subsample size for BO objective (fixed across all evals)
Nsfreq = max(1, floor(length(freq) * NsfreqFrac));

% Define deterministic objective for BO (CRN + robust replicate mean)
objectiveFunc = @(lambda) modelObjective_varred(lambda, ...
    Ne, Nr, Nw, freq, Cross, BayesIter_Conn, K, D, C, x0, Lipschitz, ...
    R_reps, repSeeds, trimFrac, Nsfreq);

% Domain in log scale (as you had)
ld_domain = optimizableVariable('l1', [lambda_space(1,1), lambda_space(1,2)], 'Transform', 'log'); % D
lc_domain = optimizableVariable('l2', [lambda_space(2,1), lambda_space(2,2)], 'Transform', 'log'); % C

% Configure BO: deterministic objective now, GP can trust observations
if indx_parallel
    results = bayesopt(objectiveFunc, [ld_domain, lc_domain], ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'MaxObjectiveEvaluations', BayesIter_Delay, ...
        'IsObjectiveDeterministic', false, ...
        'NumSeedPoints', 5, ...
        'UseParallel', true, ...
        'Verbose', 0, ...
        'PlotFcn', []);
else
    results = bayesopt(objectiveFunc, [ld_domain, lc_domain], ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'MaxObjectiveEvaluations', BayesIter_Delay, ...
        'IsObjectiveDeterministic', true, ...
        'NumSeedPoints', 5, ...
        'UseParallel', false, ...
        'Verbose', 0, ...
        'PlotFcn', []);
end

% BO optimum (according to the low-variance objective)
lambda_best = [results.XAtMinObjective.l1, results.XAtMinObjective.l2];


lambda_opt_dc =lambda_best;

end


% ======================= Objective (variance-reduced) =======================
function [Fval] = modelObjective_varred(lambda, Ne, Nr, Nw, freq, Cross, BayesIter, K, D, C, x0, Lipschitz, R_reps, repSeeds, trimFrac, Nsfreq)
% Deterministic objective via:
% - fixed frequency subset size (Nsfreq), deterministic selection
% - R_reps replicates with fixed seeds (CRN)
% - trimmed-mean aggregation of FSmooth across replicates

import functions.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.Regularization.*
import functions.auxx.TOperator.*
import functions.StochasticFISTA.*
import functions.auxx.RegSpace.*
import functions.auxx.StochasticEval.*

% Extract params from the table-like 'lambda'
lambda1 = lambda.l1;
lambda2 = lambda.l2;

% Build compact model with *scaled* D, C for this lambda
params.Dimensions.Ne = Ne;
params.Dimensions.Nr = Nr;
params.Dimensions.Nv = Nr;
params.Dimensions.Nw = Nw;
params.Data.Cross     = Cross;
params.Data.freq      = freq;
params.Compact_Model.K = K;
params.Compact_Model.D = lambda1 * D;
params.Compact_Model.C = lambda2 * C;
params.Parallel.T = 0;

% Evaluate operator T for this lambda ONCE
T = Teval(params);
% clear params; % (optional)

% Fixed, deterministic sampling across lambdas
index_stoch        = 0;                 % no random selection
index_parall_fista = 0;
Nrand              = 1;                 % single random restart per replicate

% Replicate loop with CRN seeds
FS = zeros(R_reps, 1);
for r = 1:R_reps
    % Set replicate-specific but globally fixed seed
    rng(repSeeds(r), 'twister');

    % Run inner solver (deterministic given seed + Nsfreq + index_stoch=0)
    [x_opt, ~] = stoch_fista_global([0,0,0], Ne, Nr, T, freq, ...
        index_stoch, index_parall_fista, Nsfreq, Cross, Nrand, Lipschitz, x0);

    FS(r) = x_opt.FSmooth;  % objective value to aggregate
end

% Robust aggregation (trimmed mean)
Fval = local_trimmed_mean(FS, trimFrac);
end


% Helper to call objective when we already have [l1 l2] numeric pair
function [Fval] = modelObjective_varred_tableInput(lpair, Ne, Nr, Nw, freq, Cross, BayesIter, K, D, C, x0, Lipschitz, R_reps, repSeeds, trimFrac, Nsfreq)
lambda_tbl = table(lpair(1), lpair(2), 'VariableNames', {'l1','l2'});
Fval = modelObjective_varred(lambda_tbl, Ne, Nr, Nw, freq, Cross, BayesIter, K, D, C, x0, Lipschitz, R_reps, repSeeds, trimFrac, Nsfreq);
end


% ============================ Small utility ================================
function m = local_trimmed_mean(x, trimFrac)
% Pure-MATLAB trimmed mean to avoid toolbox dependency.
% trimFrac is total fraction to trim (e.g., 0.15 trims 7.5% per tail).
x = sort(x(:));
n = numel(x);
k = floor((trimFrac * n) / 2);
if 2*k >= n
    % Degenerate case: not enough samples to trim; fall back to mean
    m = mean(x);
else
    m = mean(x(k+1 : n-k));
end
end

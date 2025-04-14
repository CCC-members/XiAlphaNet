function [best_lambda, best_index, Eval] = find_best_lambda(...
    freq, T, Cross, stoch1, stoch2, Nsfreq, x0, ...
    Ne, Nr, Nv, N_spaces, index_parall_bayes, ...
    Nrand1, Nrand2, Lipschitz, conn_delay)

import app.*
import app.functions.*
import functions.*
import functions.StochasticFISTA.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.CrossSpectrum.*
import functions.auxx.ModelVectorization.*
import functions.auxx.Regularization.*
import functions.auxx.TOperator.*
import functions.auxx.GenerateSourceSample.*
import functions.auxx.RegSpace.*
import functions.auxx.Simulations.*
import tools.*

% FIND_BEST_LAMBDA - Selects optimal lambda using Bayesian optimization and FISTA evaluation
%
% Inputs:
%   freq               - Frequency vector
%   T                  - Tensor field or transfer function structure
%   Cross              - Cross-spectrum matrix
%   stoch1             - Stochastic index set for lambda optimization
%   stoch2             - Stochastic index set for FISTA evaluation
%   Nsfreq             - Number of frequency samples
%   x0                 - Initial solution guess
%   Ne, Nr, Nv         - Number of electrodes, regions, and voxels
%   N_spaces           - Number of lambda spaces to evaluate
%   index_parall_bayes - Flag for Bayesian optimization parallelization
%   Nrand1, Nrand2     - Number of random samples for optimization and FISTA
%   Lipschitz          - Lipschitz constant
%   conn_delay         - If 1, use parfor for parallel execution
%
% Outputs:
%   best_lambda        - Optimal lambda selected
%   best_index         - Index of best lambda
%   Eval               - Array of smoothed objective values for each lambda

% Generate regularization space
lambda_space = lambda_regspace(freq, T, Cross, stoch1, Nsfreq, x0);

% Initialize evaluation containers
Eval = zeros(1, N_spaces);
Lambda = zeros(N_spaces, length(lambda_space));

% Choose loop type based on conn_delay
if conn_delay == 1
    parfor j = 1:N_spaces
        lambda_j = bayesianOptSearch(lambda_space, Ne, Nr, T, freq, stoch1, ...
            0, 0, Nsfreq, Cross, ...
            1, Lipschitz, 30, x0);
        Lambda(j,:) = lambda_j;

        [x_j, ~] = stoch_fista_global(lambda_j, Ne, Nv, T, freq, ...
            stoch2, 0, Nsfreq, Cross, ...
            1, Lipschitz, x0);

        Eval(j) = x_j.FSmooth;
    end
else
    for j = 1:N_spaces
        lambda_j = bayesianOptSearch(lambda_space, Ne, Nr, T, freq, stoch1, ...
            0, 0, Nsfreq, Cross, ...
            1, Lipschitz, 1, x0);
        Lambda(j,:) = lambda_j;

        [x_j, ~] = stoch_fista_global(lambda_j, Ne, Nv, T, freq, ...
            stoch2, 0, Nsfreq, Cross, ...
            30, Lipschitz, x0);

        Eval(j) = x_j.FSmooth;
    end
end

% Identify best lambda
[~, best_index] = min(Eval);
best_lambda = Lambda(best_index, :);
end

function [lambda_opt] = bayesianOptSearch(lambda_space, Ne, Nv, T, freq, ...
    index_stoch, index_parall_fist, index_parall_bayes, Nsfreq, Cross, ...
    Nrand, Lipschitz, BayesIter, x0)
% BAYESIANOPTSEARCH
% Uses Bayesian Optimization to find optimal regularization parameters
% for an inverse solution model with Lasso constraints.
%
% Inputs:
%   lambda_space       - Vector [max_l1, max_l2, max_l3] specifying the 
%                        maximum range for lambda parameters.
%   Ne, Nv, T, freq    - Model dimensions and frequency information.
%   index_stoch        - Flag for stochastic evaluation of the objective.
%   index_parall_fist  - Flag for parallelization inside FISTA solver.
%   index_parall_bayes - Flag for parallelization inside bayesopt.
%   Nsfreq             - Number of spectral frequencies.
%   Cross              - Cross-spectral data structure.
%   Nrand              - Number of random seeds / replications.
%   Lipschitz          - Lipschitz constant for solver.
%   BayesIter          - Maximum number of Bayesian optimization iterations.
%   x0                 - Initial guess for solver.
%
% Outputs:
%   lambda_opt - Optimal lambda values (vector of length 3).
%

    %% Define the objective function to be minimized
    objectiveFunc = @(lambda) modelObjective(lambda, Ne, Nv, T, freq, ...
                                             index_stoch, index_parall_fist, ...
                                             Nsfreq, Cross, Nrand, ...
                                             Lipschitz, x0);

    %% Define Bayesian optimization search space
    l1_domain = optimizableVariable('l1', [1e-10, lambda_space(1)], 'Transform', 'log');
    l2_domain = optimizableVariable('l2', [1e-10, lambda_space(2)], 'Transform', 'log');
    l3_domain = optimizableVariable('l3', [1e-10, lambda_space(3)], 'Transform', 'log');

    %% Run Bayesian Optimization
    if index_parall_bayes == 1
        results = bayesopt(objectiveFunc, [l1_domain, l2_domain, l3_domain], ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxObjectiveEvaluations', BayesIter, ...
            'IsObjectiveDeterministic', false, ...
            'NumSeedPoints', 12, ...
            'Verbose', 0, ...
            'PlotFcn', [], ...
            'UseParallel', true);
    else
        results = bayesopt(objectiveFunc, [l1_domain, l2_domain, l3_domain], ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxObjectiveEvaluations', BayesIter, ...
            'IsObjectiveDeterministic', false, ...
            'NumSeedPoints', 12, ...
            'Verbose', 0, ...
            'PlotFcn', []);
    end

    %% Extract optimal lambda values
    lambda_opt = [results.XAtMinObjective.l1, ...
                  results.XAtMinObjective.l2, ...
                  results.XAtMinObjective.l3];

end


%% Nested model objective function
function [bic_val, solution] = modelObjective(lambda, Ne, Nv, T, freq, ...
    index_stoch, index_parall, Nsfreq, Cross, Nrand, Lipschitz, x0)
% MODELOBJECTIVE
% Wrapper that computes the BIC objective for given lambda values.
%
% Inputs/Outputs as in bayesianOptSearch.

    import functions.*
    import functions.auxx.*
    import functions.auxx.BayesOptimization.*

    % Compute BIC
    BIC_result = BIC(lambda, Ne, Nv, T, freq, ...
                     index_stoch, index_parall, ...
                     Nsfreq, Cross, Nrand, Lipschitz, x0);

    % Objective value to minimize
    bic_val  = BIC_result.Feval;
    solution = BIC_result.Solution;
end

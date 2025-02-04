function [lambda_opt_dc] = bayes_search_conn_delay(lambda_space, parameters)
    
    % Define the objective function to be minimized using Bayesian Optimization
    objectiveFunc = @(lambda) modelObjective(lambda, parameters);
    
    % Set up the domain for lambda parameters
    ld_domain = optimizableVariable('l1', [lambda_space(1,1), lambda_space(1,2)], 'Transform', 'log'); % s
    lc_domain = optimizableVariable('l2', [lambda_space(2,1), lambda_space(2,2)], 'Transform', 'log'); % e
   
    
    % Run Bayesian Optimization
    tic;
    if parameters.Parallel.Conn_Delay == 1
        results = bayesopt(objectiveFunc, [ld_domain,lc_domain], ...
                           'AcquisitionFunctionName', 'expected-improvement-plus', ...
                           'MaxObjectiveEvaluations', parameters.BayesIter0, ...
                           'IsObjectiveDeterministic', false, ...
                           'NumSeedPoints', 5, ...
                           'Verbose', 0,'UseParallel',true);
    else 
       results = bayesopt(objectiveFunc, [ld_domain,lc_domain], ...
                       'AcquisitionFunctionName', 'expected-improvement-plus', ...
                       'MaxObjectiveEvaluations', parameters.BayesIter0, ...
                       'IsObjectiveDeterministic', false, ...
                       'NumSeedPoints', 5, ...
                       'Verbose', 0,'UseParallel',false);
    end
    toc;
    
    % Extract optimal lambda values and the corresponding solution
    lambda_opt_dc(1) = results.XAtMinObjective.l1;
    lambda_opt_dc(2) = results.XAtMinObjective.l2;

    close Figure 1
    close Figure 2
end

% Model objective function that computes the AIC and solution for given lambda
function [LL_val] = modelObjective(lambda, parameters)
    % Extract lambda values from the table using variable names
    lambda1 = lambda.l1;
    lambda2 = lambda.l2;
    
    % Use the lambda values in your calculations
    temp_parameters = parameters; % Create a copy to avoid modifying the original
    temp_parameters.Compact_Model.D = lambda1 * temp_parameters.Compact_Model.D;
    temp_parameters.Compact_Model.C = lambda2 * temp_parameters.Compact_Model.C;
    temp_parameters = Teval(temp_parameters);
    
    % Estimate the number of frequencies that approximate the f and dF with a relative error less than 5%
    %k_min = findMinimumK(temp_parameters, 10, 5);
    temp_parameters.Stochastic.stoch = 1;
    temp_parameters.Stochastic.Nsfreq =25;
    temp_parameters.Stochatic.Niter = 1;
    temp_parameters = sample_frequencies(temp_parameters);
    temp_parameters.Threshold = activation_threshold(temp_parameters);
    temp_parameters.Lipschitz = 0.01; % Estimate or assign the Lipschitz constant
    lambda_space = [100, 1000, 1000]; % Adjust as needed
    disp('-->> Initializing Stochastic FISTA global optimizer...')
    [lambda_opt] = bayesianOptSearch(lambda_space, temp_parameters);
    temp_parameters.Stochastic.stoch = 0;
    temp_parameters.Stochastic.Niter = 1;
    temp_parameters = sample_frequencies(temp_parameters);
    tic
    [x_opt, ~] = stoch_fista_global(lambda_opt, temp_parameters);
    toc
    LL_val = x_opt.FSmooth;
end

function [lambda_opt_dc] = bayes_search_conn_delay(lambda_space, Ne,Nr,Nw,freq,Cross,BayesIter_Conn,K,D,C,indx_parallel,BayesIter_Delay)
    
    % Define the objective function to be minimized using Bayesian Optimization
    objectiveFunc = @(lambda) modelObjective(lambda,Ne,Nr,Nw,freq,Cross,BayesIter_Conn,K,D,C);
    
    % Set up the domain for lambda parameters
    ld_domain = optimizableVariable('l1', [lambda_space(1,1), lambda_space(1,2)], 'Transform', 'log'); % s
    lc_domain = optimizableVariable('l2', [lambda_space(2,1), lambda_space(2,2)], 'Transform', 'log'); % e
   
    
    % Run Bayesian Optimization
    tic;
    if indx_parallel
        results = bayesopt(objectiveFunc, [ld_domain,lc_domain], ...
                           'AcquisitionFunctionName', 'expected-improvement-plus', ...
                           'MaxObjectiveEvaluations', BayesIter_Delay, ...
                           'IsObjectiveDeterministic', false, ...
                           'NumSeedPoints', 5, ...
                           'Verbose', 0,'UseParallel',true);
    else 
       results = bayesopt(objectiveFunc, [ld_domain,lc_domain], ...
                       'AcquisitionFunctionName', 'expected-improvement-plus', ...
                       'MaxObjectiveEvaluations', BayesIter_Delay, ...
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
function [LL_val] = modelObjective(lambda,Ne,Nr,Nw,freq,Cross,BayesIter,K,D,C)
    % Extract lambda values from the table using variable names
    lambda1 = lambda.l1;
    lambda2 = lambda.l2;
    temp_parameters.Dimensions.Ne = Ne;
    temp_parameters.Dimensions.Nr=Nr;
    temp_parameters.Dimensions.Nv=Nr;
    temp_parameters.Dimensions.Nw=Nw;
    temp_parameters.Data.Cross = Cross;
    temp_parameters.Data.freq = freq;
    temp_parameters.Compact_Model.K = K;
    temp_parameters.Compact_Model.D = D;
    temp_parameters.Compact_Model.C = C;
    % Use the lambda values in your calculations
    %temp_parameters = parameters; % Create a copy to avoid modifying the original
    temp_parameters.Compact_Model.D = lambda1 * temp_parameters.Compact_Model.D;
    temp_parameters.Compact_Model.C = lambda2 * temp_parameters.Compact_Model.C;
    temp_parameters.Parallel.T=0;
    T = Teval(temp_parameters);
    clear temp_parameters;
    % Estimate the number of frequencies that approximate the f and dF with a relative error less than 5%
    k_min = 30;%findMinimumK(freq,T,Cross, 5, 20,0);
    index_stoch = 1;
    index_parall_fist = 0;
    index_parall_bayes= 0;
    Nsfreq = k_min;
    Nrand = 1;
    Lipschitz = 0.01;
    lambda_space = [100, 1000, 1000]; % Adjust as needed
    disp('-->> Initializing Stochastic FISTA global optimizer...')
    [lambda_opt] = bayesianOptSearch(lambda_space,Ne,Nr,T,freq,index_stoch,index_parall_fist,index_parall_bayes,Nsfreq,Cross,Nrand,Lipschitz,BayesIter);
    index_stoch = 0;
    Nrand = 1;
    tic
    [x_opt, ~] = stoch_fista_global(lambda_opt, Ne,Nr,T,freq,index_stoch,index_parall_fist,Nsfreq,Cross,Nrand,Lipschitz);
    toc
    LL_val = x_opt.FSmooth;
end

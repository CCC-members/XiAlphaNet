function [lambda_opt] = bayesianOptSearch(lambda_space,Ne,Nv,T,freq,index_stoch,index_parall_fist,index_parall_bayes,Nsfreq,Cross,Nrand,Lipschitz,BayesIter)
% BAYESIANOPTSEARCH Uses Bayesian Optimization to find optimal regularization
% parameters for an inverse solution model with Lasso constraints.
%
% Inputs:
%   lambda_space - Struct with fields 'l1' and 'l2' specifying the range for lambda parameters.
%   parameters   - Struct containing other parameters for the model.
%
% Outputs:
%   lambda_opt   - Optimal lambda values found (struct with fields 'l1' and 'l2').
%   x_opt        - Optimal solution corresponding to the optimal lambda.

%% Stocastic Exploration of the Regularization space 
% Define the objective function to be minimized using Bayesian Optimization
    objectiveFunc = @(lambda) modelObjective(lambda, Ne,Nv,T,freq,index_stoch,index_parall_fist,Nsfreq,Cross,Nrand,Lipschitz);
    
    % Set up the domain for lambda parameters
    l1_domain = optimizableVariable('l1', [10^(-10), lambda_space(1)], 'Transform', 'log'); % s
    l2_domain = optimizableVariable('l2', [20, lambda_space(2)], 'Transform', 'log'); % e
    l3_domain = optimizableVariable('l3', [20, lambda_space(3)], 'Transform', 'log'); % a
    
    % Run Bayesian Optimization
    if index_parall_bayes == 1
       results = bayesopt(objectiveFunc, [l1_domain, l2_domain, l3_domain], ...
                       'AcquisitionFunctionName', 'expected-improvement-plus', ...
                       'MaxObjectiveEvaluations', BayesIter, ...
                       'IsObjectiveDeterministic', false, ...
                       'NumSeedPoints', 5, ...
                       'Verbose', 0,'PlotFcn',[],'UseParallel',true);

    else 
        results = bayesopt(objectiveFunc, [l1_domain, l2_domain, l3_domain], ...
                       'AcquisitionFunctionName', 'expected-improvement-plus', ...
                       'MaxObjectiveEvaluations', BayesIter, ...
                       'IsObjectiveDeterministic', false, ...
                       'NumSeedPoints', 5, ...
                       'Verbose', 0,'PlotFcn',[]);
    end 

    % Extract optimal lambda values and the corresponding solution
    lambda_opt(1) = results.XAtMinObjective.l1;
    lambda_opt(2) = results.XAtMinObjective.l2;
    lambda_opt(3) = results.XAtMinObjective.l3;

%  close Figure 2
% % close Figure 2
end
% Model objective function that computes the AIC and solution for given lambda
function [bic_val, solution] = modelObjective(lambda, Ne,Nv,T,freq,index_stoch,index_parall,Nsfreq,Cross,Nrand,Lipschitz)
import functions.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
% Assuming AIC function exists and returns a structure with fields .Feval and .Solution
BIC_result = BIC(lambda, Ne,Nv,T,freq,index_stoch,index_parall,Nsfreq,Cross,Nrand,Lipschitz);
bic_val = BIC_result.Feval; % Minimize this value
solution = BIC_result.Solution; % Save solution for output
end

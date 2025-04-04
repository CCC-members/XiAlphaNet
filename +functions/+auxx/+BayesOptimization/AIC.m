function [AICV] = AIC(lambda, parameters)
% AIC Computes the Akaike Information Criterion for a given model configuration.
%
% Inputs:
%   lambda - A vector  specifying regularization strengths
%            for alpha and xi components, respectively.
%   parameters - A struct containing various settings and hyperparameters for the model,
%                including the connectivity matrix C, dimension of ROI, and the number of 
%                frequency bands (Nsfreq).
%
% Outputs:
%   AICV - A struct with fields 'Solution', which holds the optimized model parameters,
%          and 'Feval', which is the computed AIC value for the model.

% Set regularization parameters based on input

% Execute FISTA to find the optimal parameters
lambda = table2array(lambda);
[X_opt,~] = stoch_fista_global(lambda,parameters);

% Calculate the degrees of freedom for regularization
dgf = nnz(X_opt.Solution);
f = X_opt.FSmooth;
% Calculate AIC value
AICV.Solution = X_opt;     % Save optimized parameters
AICV.Feval = 2*dgf + 2*f ;  % Sum degrees of freedom and model's objective function value as AIC
end

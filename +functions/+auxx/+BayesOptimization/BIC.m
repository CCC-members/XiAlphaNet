function [BICV] = BIC(lambda, Ne,Nv,T,freq,index_stoch,index_parll,Nsfreq,Cross,Nrand,Lipschitz)
import functions.*
import functions.auxx.*
import functions.StochasticFISTA.*
import functions.auxx.ModelVectorization.*
% BIC Computes the Bayesian Information Criterion for a given model configuration.
%
% Inputs:
%   lambda - A vector  specifying regularization strengths
%            for alpha and xi components, respectively.
%   parameters - A struct containing various settings and hyperparameters for the model,
%                including the connectivity matrix C, dimension of ROI, and the number of 
%                frequency bands (Nsfreq).
%
% Outputs:
%   BICV - A struct with fields 'Solution', which holds the optimized model parameters,
%          and 'Feval', which is the computed AIC value for the model.

% Set regularization parameters based on input

% Execute FISTA to find the optimal parameters
lambda = table2array(lambda);
[X_opt,~] = stoch_fista_global(lambda, Ne,Nv,T,freq,index_stoch,index_parll,Nsfreq,Cross,Nrand,Lipschitz);
N  = length(freq);
[e,a,s2] = x2v(X_opt.Solution);
epsilon=1/length(X_opt.Solution);
% Calculate the degrees of freedom for regularization
dgf = nnz(X_opt.Solution);%+1/(norm(a,"fro")+epsilon)+1/(norm(e,"fro")+epsilon);
f = X_opt.FSmooth;
% Calculate BIC value
BICV.Solution = X_opt;     % Save optimized parameters
BICV.Feval = log(N)*dgf + 2*f ;  % Sum degrees of freedom and model's objective function value as AIC
end

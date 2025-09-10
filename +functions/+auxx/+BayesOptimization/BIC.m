function [BICV] = BIC(lambda, Ne, Nv, T, freq, index_stoch, index_parll, Nsfreq, Cross, Nrand, Lipschitz, x0)
import functions.*
import functions.auxx.*
import functions.StochasticFISTA.*
import functions.auxx.ModelVectorization.*
% BIC Computes the Bayesian Information Criterion for a given model configuration.
%
% Inputs:
%   lambda - Regularization vector (can be numeric, table row from bayesopt, or struct with fields l1,l2,l3)
%   Ne,Nv,T,freq,index_stoch,index_parll,Nsfreq,Cross,Nrand,Lipschitz,x0 - pipeline parameters
%
% Outputs:
%   BICV - struct with fields:
%          .Solution : optimized model output from stoch_fista_global
%          .Feval    : scalar BIC value to be minimized

% --- Coerce lambda to numeric row vector (handles table/struct inputs)
if istable(lambda)
    lambda = table2array(lambda(1,:));
elseif isstruct(lambda) && all(isfield(lambda, {'l1','l2','l3'}))
    lambda = [lambda.l1, lambda.l2, lambda.l3];
end
lambda = double(lambda(:)).';

% --- Solve the model at this lambda
[X_opt, ~] = stoch_fista_global(lambda, Ne, Nv, T, freq, ...
    index_stoch, index_parll, Nsfreq, Cross, Nrand, Lipschitz, x0);

% --- Data-fit term (must be the pure data term, no regularization penalties)
f_data = X_opt.FSmooth;

% --- Degrees of freedom (effective nonzeros; exclude nuisance if split available)
sol = X_opt.Solution;
tol = 1e-8;
try
    [e, a, ~] = x2v(sol);                             % exclude s2 (noise variance) if present
    df = nnz(abs(e) > tol) + nnz(abs(a) > tol);
catch
    df = nnz(abs(sol(:)) > tol);
end

% --- Effective sample size n (match your data-fit construction)
Nf = numel(freq);
n_eff = Nf * Ne;                                     % conservative fallback (e.g., diagonal/magnitude per freq)
if ~isempty(Cross)
    sz = size(Cross);
    if numel(sz) >= 2 && sz(1) == Ne && sz(2) == Ne  % full (Hermitian) cross-spectra available
        n_eff = Nf * (Ne * (Ne + 1) / 2);           % unique real dof per frequency
    end
end
n_eff = max(1, n_eff);

% --- BIC value
BICV.Solution = X_opt;
BICV.Feval    = log(n_eff) * df + 2 * f_data;
end

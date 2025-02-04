function X = prox_ridge(U, opts)
% function X = prox_ridge(U, opts)
% Description:
%   Proximal operator for the ridge (L2) regularization. Solves the following problem:
%       U is a matrix of size d x k 
%       lambda is a positive scalar
%       X = argmin_X 0.5*||X - U||_F^2 + lambda * ||X||_F^2 
%       if `opts.pos = true`, then this function solves the problem with 
%       positive constraints.
% Inputs: 
%   U: double dense matrix d x k 
%   opts: struct 
%       opts.lambda: regularization parameter (scalar)
%       opts.pos: positive constraint (default = false)
% Outputs: 
%   X: a full matrix of size d x k
% -----------------------------------------------
% Author: Adapted from Tiep Vu's code
% -----------------------------------------------
   
    lambda = opts;

    % Compute the proximal operator for ridge regularization
    X = U / (1 + 2 * lambda);
    X = max(0, X);
   
end

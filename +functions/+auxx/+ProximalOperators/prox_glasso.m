function X = prox_glasso(U, opts)
% function X = proj_group_lasso(U, opts)
% Description:
%   Proximal operator for the Group Row Lasso. Solves the following problem:
%       U is a matrix of size d x k 
%       lambda is a positive scalar or a column vector
%       X = argmin_X 0.5*||X - U||_F^2 + lambda * sum_i ||X_i||_2 
%       where X_i and U_i are the i-th rows of X and U respectively.
%       if `opts.pos = true`, then this function solves the problem with 
%       positive constraints.
% Inputs: 
%   U: double dense matrix d x k 
%   opts: struct 
%       opts.lambda: regularization parameter (scalar or column vector)
%       opts.pos: positive constraint (default = false)
% Outputs: 
%   X: a full matrix of size d x k
% -----------------------------------------------
% Author: Adapted from Tiep Vu's code
% -----------------------------------------------
    [N,~]  = size(U);
    lambda = opts/N;

    % Compute the L2 norm of each row of U
    row_norms = sqrt(sum(U.^2, 2));
    
    % Compute the scaling factors
    scaling_factors = max(0, 1 - lambda ./ row_norms);
    
    % Ensure scaling_factors is a column vector
    scaling_factors = scaling_factors(:);
    
    % Apply the scaling factors to each row
    X = U .* scaling_factors;
   
    X = max(0, X);

end

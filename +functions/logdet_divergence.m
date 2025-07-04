function D = logdet_divergence_hermitian(A, B)
%LOGDET_DIVERGENCE_HERMITIAN Computes the LogDet divergence for Hermitian matrices
%
%   D = logdet_divergence_hermitian(A, B)
%
%   Inputs:
%       A, B - Hermitian positive definite or semi-definite matrices
%   Output:
%       D - LogDet divergence

    % Check size
    if ~isequal(size(A), size(B))
        error('Input matrices must be the same size.');
    end
  A = 0.5*(A+A');
  B = 0.5*(B+B');
    % Check Hermitian
    if ~ishermitian(A) || ~ishermitian(B)
        error('Both matrices must be Hermitian.');
    end

    % Check if B is invertible
    % if rcond(B) < 1e-12
    %     error('Matrix B is nearly singular or not positive definite.');
    % end

    % Compute B^{-1}A
    BinvA = pinv(B)* A+0.0001*eye(size(A));

    % Use eigenvalue method for log-det to handle complex/Hermitian
    lambda = eig(BinvA);
    
    % Ensure eigenvalues are positive (required for log)
    if any(real(lambda) <= 0)
        warning('Matrix B^{-1}A has non-positive eigenvalues. Divergence may be complex or undefined.');
    end

    % Compute trace and log-det via eigenvalues
    D = real(trace(BinvA) - sum(log(lambda)) - length(lambda));
end

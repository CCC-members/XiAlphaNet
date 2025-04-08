function S = generate_complex_wishart(Sigma, n)
% Generate a simulated covariance matrix from a complex Wishart distribution.
%
% Parameters:
%   Sigma - The p x p complex covariance matrix (must be positive definite).
%   n     - Degrees of freedom (number of samples).
%
% Returns:
%   S     - A simulated p x p covariance matrix from the complex Wishart distribution.
import functions.*
import functions.auxx.*
import functions.auxx.DataPreprosessing.*
import functions.auxx.OptimizedOperations.*
% Ensure that Sigma is a square matrix
Sigma  = regularize_tensor(Sigma);
[p, p_col] = size(Sigma);
if p ~= p_col
    error('Sigma must be a square matrix.');
end

% Cholesky decomposition of Sigma (lower triangular matrix)
L = chol(Sigma, 'lower');

% Generate a p x n matrix of standard complex normal random variables
% Each entry is of the form (1/sqrt(2)) * (a + i*b),
% where a and b are independent N(0,1) real random variables.
W_real = randn(p, n) / sqrt(2);
W_imag = randn(p, n) / sqrt(2);
W = W_real + 1i * W_imag;

% Compute Z = L * W
Z = L * W;

% Compute the sample covariance matrix S = Z * Z^H
% Note: In MATLAB, the apostrophe operator (') performs the conjugate transpose for complex matrices
S = Z * Z';
end

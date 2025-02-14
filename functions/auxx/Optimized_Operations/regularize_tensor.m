function [A_reg] = regularize_tensor(A)
% This function regularizes each slice of the tensor A using Tikhonov regularization.
% It finds the optimal lambda for each slice to make it positive definite.
% If the slice is already positive definite, no regularization is applied.
%
% Input:
%   A - 3D tensor of size (m, n, p), where each A(:,:,j) is a matrix.
%       If A is 2D, it will be treated as a single slice.
%
% Output:
%   A_reg - The regularized 3D tensor (or 2D matrix if A is 2D).
%   lambda_final - The maximum lambda used for regularization across all slices.

% Check if A is a 2D matrix or 3D tensor
if ndims(A) == 2
    % If A is 2D, treat it as a single slice
    A = reshape(A, size(A, 1), size(A, 2), 1);
end

% Get the size of the input tensor
[m, n, p] = size(A);

% Initialize the vector to store the optimal lambdas for each slice
lambda_optimal = zeros(1, p);

% Regularization parameter range or starting value
lambda_initial = 1e-6;  % Start with a small lambda value
lambda_increase_factor = 1.5;  % Factor to increase lambda if needed
max_iterations = 50;  % Max iterations to search for positive definiteness
epsilon = 1e-6;  % Small value to add to the matrix to ensure positive definiteness

% Loop through each slice A(:,:,j)
for j = 1:p
    A_j = A(:,:,j); % Extract the j-th slice

    % Ensure the matrix is symmetric if needed
    A_j = (A_j + A_j') / 2;

    % Check if the matrix is already positive definite
    eigvals = eig(A_j);
    min_eigenvalue = min(eigvals);

    if min_eigenvalue > 0
        % If the matrix is already positive definite, no regularization is needed
        lambda_optimal(j) = 0;  % No regularization applied
    else
        % Otherwise, initialize lambda for this slice and the iteration counter
        lambda = lambda_initial;
        iteration = 0;

        % Loop to find the lambda that makes the matrix positive definite
        while iteration < max_iterations
            % Apply Tikhonov regularization: A_j + lambda * I
            A_j_reg = A_j + lambda * eye(m); % Regularization term

            % Compute the smallest eigenvalue of the regularized matrix
            eigvals_reg = eig(A_j_reg);
            min_eigenvalue_reg = min(eigvals_reg);

            % Check if the matrix is positive definite
            if min_eigenvalue_reg > 0
                % If positive definite, stop the loop and record lambda
                lambda_optimal(j) = lambda;
                break;
            else
                % Otherwise, increase lambda and continue the search
                lambda = lambda * lambda_increase_factor;
                iteration = iteration + 1;
            end
        end

        % If the matrix is still not positive definite after max iterations
        if iteration == max_iterations
            disp(['Warning: Slice ', num2str(j), ' could not be made positive definite']);
            % Find the minimum eigenvalue and make it positive by adding epsilon
            min_eigenvalue_reg = min(eigvals);
            A_j = A_j + abs(min_eigenvalue_reg) + epsilon;  % Make it positive definite
            lambda_optimal(j) = abs(min_eigenvalue_reg) + epsilon;  % Record the regularization term
        end
    end
end

% The final regularization parameter is the maximum of the optimal lambdas
lambda_final = max(lambda_optimal);

% Apply the final regularization to all slices with the same lambda
for j = 1:p
    A(:,:,j) = A(:,:,j) + lambda_final * eye(m); % Regularize each slice with lambda_final
end

% If A was originally a 2D matrix, return the regularized 2D matrix
A_reg = A;
if ndims(A_reg) == 3 && size(A_reg, 3) == 1
    A_reg = A_reg(:,:,1);  % Convert back to 2D if it was originally 2D
end
end

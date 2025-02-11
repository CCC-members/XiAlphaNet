function log_det_A = logdet(A, epsilon)
    % Check if epsilon is provided, otherwise set a default small value
    if nargin < 2
        epsilon = 1e-6;  % Default regularization term
    end

    % Add regularization term to diagonal
    A_reg = A + epsilon * eye(size(A));

    % Compute the eigenvalues of the regularized matrix
    eigenvalues = eig(A_reg);

    % Calculate the logarithm of the determinant
    % Remove non-positive eigenvalues before taking log to avoid complex results
    positive_eigenvalues = eigenvalues(eigenvalues > 0);
    log_det_A = sum(log(positive_eigenvalues));
   
end

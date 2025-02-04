function lambda = compute_lambda(gradient_matrix, L)
    % Function to compute the regularization parameter lambda
    % given the gradient matrix and the Lipschitz constant L.
    %
    % Parameters:
    % gradient_matrix: A matrix where each row represents a group gradient
    % L: The Lipschitz constant
    %
    % Returns:
    % lambda: The computed regularization parameter

    % Compute the L2 norm of each row
    row_norms = sqrt(sum(gradient_matrix.^2, 2));
    
    % Find the maximum L2 norm
    max_row_norm = max(row_norms);
    
    % Compute lambda
    lambda = (1/L) * max_row_norm;
end

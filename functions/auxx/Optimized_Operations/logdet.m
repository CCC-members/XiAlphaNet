function log_det_A = logdet(A, epsilon)
    % Check if epsilon is provided, otherwise set a default small value
    if nargin < 2
        epsilon = 1e-6;  % Default regularization term
    end

    % Check the condition number of the matrix to decide whether to regularize more
    cond_A = cond(A);
    if cond_A > 1e12  % If the matrix is ill-conditioned, increase the regularization
        epsilon = epsilon * 10;
    end

    % Regularize the matrix (use a larger epsilon for ill-conditioned matrices)
    A_reg = A + epsilon * eye(size(A));

    % Check if the matrix is symmetric and positive semi-definite
    while true
        % Check the eigenvalues to see if matrix is positive definite
        eigvals = eig(A_reg);
        
        % If all eigenvalues are positive, proceed
        if all(eigvals > 0)
            break;
        else
            % Increase epsilon until the matrix becomes positive definite
            epsilon = epsilon * 2;  % Double epsilon each time
            A_reg = A + epsilon * eye(size(A));  % Re-regularize the matrix
            warning('Matrix not positive definite, increasing epsilon to %.5f', epsilon);
        end
    end

    % Check if the matrix is symmetric and positive semi-definite
    if isequal(A_reg, A_reg') && all(eig(A_reg) > 0)
        % Use Cholesky Decomposition for symmetric positive semi-definite matrices
        try
            R = chol(A_reg);  % Cholesky decomposition
            log_det_A = 2 * sum(log(diag(R)));  % log(det(A)) = 2 * sum(log(diag(R)))
        catch
            % If Cholesky fails, fall back to SVD
            warning('Cholesky failed, using SVD for log-determinant.');
            [U, S, V] = svd(A_reg);
            log_det_A = sum(log(diag(S)));  % Use the singular values from SVD
        end
    else
        % Use SVD for general matrices
        [U, S, V] = svd(A_reg);
        log_det_A = sum(log(diag(S)));  % Use the singular values from SVD
    end
end

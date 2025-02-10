function X = prox_group_lasso(U, lambda)
    % Compute the L2 norm of each row of U
    row_norms = sqrt(sum(U.^2, 2));
    
    % Compute the scaling factors
    scaling_factors = max(0, 1 - lambda ./ row_norms);
    
    % Apply the scaling factors to each row
    X = U .* scaling_factors;
end
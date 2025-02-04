function numcheck_xialpha_derivatives(K, Nr, omega)
    % Number of rows in e and a
    % Nr: number of rows in e and a
    % omega: scalar or vector
    % K: number of random trials

    % Initialize accumulators for the relative absolute errors
    xi_error_sum = 0;
    alpha_error_sum = 0;

    % Define small perturbation
    h = 1e-3;

    for k = 1:K
        % Generate random matrices e and a
        e = rand(Nr, 3);
        a = rand(Nr, 4);

        % Analytical derivatives
        [dXi_de_analytical, dAlpha_da_analytical] = dXiAlpha_deda(e, a, omega);

        % Initialize numerical derivative matrices
        dXi_de_numerical = zeros(size(e));
        dAlpha_da_numerical = zeros(size(a));

        % Numerical derivatives for xi_omega
        for i = 1:Nr
            for j = 1:3
                e_perturb = e;
                e_perturb(i, j) = e_perturb(i, j) + h;
                xi_omega_perturb = e_perturb(:,1) ./ (1 + e_perturb(:,2) .* omega.^2).^e_perturb(:,3);
                xi_omega_original = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
                dXi_de_numerical(i, j) = (xi_omega_perturb(i) - xi_omega_original(i)) / h;
            end
        end

        % Numerical derivatives for alpha_omega
        for i = 1:Nr
            for j = 1:4
                a_perturb = a;
                a_perturb(i, j) = a_perturb(i, j) + h;
                alpha_omega_perturb = a_perturb(:,1) ./ (1 + a_perturb(:,2) .* (omega - a_perturb(:,4)).^2).^a_perturb(:,3);
                alpha_omega_original = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
                dAlpha_da_numerical(i, j) = (alpha_omega_perturb(i) - alpha_omega_original(i)) / h;
            end
        end

        % Compute relative absolute error
        xi_error = abs(dXi_de_analytical - dXi_de_numerical) ./ abs(dXi_de_analytical);
        alpha_error = abs(dAlpha_da_analytical - dAlpha_da_numerical) ./ abs(dAlpha_da_analytical);

        % Handle potential NaN values due to division by zero
        xi_error(isnan(xi_error)) = 0;
        alpha_error(isnan(alpha_error)) = 0;

        % Accumulate errors
        xi_error_sum = xi_error_sum + sum(xi_error(:)) / numel(xi_error);
        alpha_error_sum = alpha_error_sum + sum(alpha_error(:)) / numel(alpha_error);
    end

    % Compute mean relative absolute error
    mean_xi_error = xi_error_sum / K;
    mean_alpha_error = alpha_error_sum / K;

    % Display results
    fprintf('Mean relative absolute error for xi_omega: %e\n', mean_xi_error);
    fprintf('Mean relative absolute error for alpha_omega: %e\n', mean_alpha_error);
end

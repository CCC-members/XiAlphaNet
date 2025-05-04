function result = fitZIG_random(Y_roi, age, max_iter, tol, verbose)
% FITZIG_FIXED - EM estimation of ZIG with only fixed effects
% OUTPUT:
%   result: struct with beta, gamma, sigma, p-values

    if nargin < 3, max_iter = 100; end
    if nargin < 4, tol = 1e-6; end
    if nargin < 5, verbose = true; end

    [Nr, N] = size(Y_roi);
    age = age(:);
    Age_long = repmat(age, Nr, 1);
    Age2_long = Age_long.^2;

    % Fixed effects design matrices
    X_mu = [ones(Nr*N, 1), Age_long, Age2_long];
    X_pi = [ones(Nr*N, 1), Age_long];

    % Flatten data
    Y_long = reshape(Y_roi', [], 1);
    is_pos = Y_long > 0;

    % Initialize parameters
    beta = X_mu(is_pos, :) \ Y_long(is_pos);
    gamma = glmfit(Age_long, ~is_pos, 'binomial', 'link', 'logit');
    sigma = std(Y_long(is_pos));

    prev_ll = -Inf;

    for iter = 1:max_iter
        % === E-step ===
        mu = X_mu * beta;
        pi_logit = X_pi * gamma;
        pi = 1 ./ (1 + exp(-pi_logit));
        pi = min(max(pi, 1e-6), 1 - 1e-6);

        phi0 = normpdf(0, mu, sigma);
        tau = zeros(Nr*N, 1);
        tau(is_pos) = 1;

        nonpos = ~is_pos;
        denom = pi(nonpos) + (1 - pi(nonpos)) .* phi0(nonpos) + eps;
        tau(nonpos) = ((1 - pi(nonpos)) .* phi0(nonpos)) ./ denom;
        tau = min(max(tau, 1e-6), 1 - 1e-6);

        % === M-step ===
        % Update fixed effects for mean
        % Now we update W diagonally using tau directly, no need to form the full W matrix
        W_diag = tau;  % This is the diagonal of W, a vector of length Nr*N

        % Update beta (weighted least squares)
        beta = (X_mu' * (W_diag .* X_mu)) \ (X_mu' * (W_diag .* Y_long));

        % Update residual variance
        residuals = Y_long - X_mu * beta;
        sigma2 = sum(W_diag .* residuals.^2) / sum(W_diag);
        sigma = sqrt(sigma2);

        % Update zero-inflation parameters
        [gamma, dev, stats] = glmfit(Age_long, 1 - tau, 'binomial', 'link', 'logit');
        se_gamma = stats.se;

        % Log-likelihood
        logL_zero = log(pi + (1 - pi) .* phi0 + eps);
        logL_pos  = log(1 - pi + eps) + log(normpdf(Y_long, mu, sigma) + eps);
        ll = sum((1 - tau) .* logL_zero) + sum(tau .* logL_pos);

        if verbose
            fprintf('Iter %3d: loglik = %.6f, mean(1-tau) = %.4f\n', iter, ll, mean(1 - tau));
        end
        if abs(ll - prev_ll) < tol, break; end
        prev_ll = ll;
    end

    % Standard errors and p-values
    XtWX = X_mu' * (W_diag .* X_mu);
    se_beta = sqrt(diag(pinv(XtWX) * sigma2));
    z_beta = beta ./ se_beta;
    p_beta = 2 * (1 - normcdf(abs(z_beta)));

    if all(~isnan(se_gamma))
        z_gamma = gamma ./ se_gamma;
        p_gamma = 2 * (1 - normcdf(abs(z_gamma)));
    else
        p_gamma = NaN(size(gamma));
    end
    mu_est = beta(1) + beta(2)*age + beta(3)*age.^2;
    pi_est = 1 ./ (1 + exp(-(gamma(1) + gamma(2)*age)));

    % === Output ===
    result.beta = beta;
    result.gamma = gamma;
    result.sigma = sigma;
    result.loglik = ll;
    result.iter = iter;
    result.p_beta = p_beta;
    result.p_gamma = p_gamma;
    result.mu_est = mu_est;
    result.pi_est = pi_est;
end
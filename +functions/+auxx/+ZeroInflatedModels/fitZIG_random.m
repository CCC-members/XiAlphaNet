function result = fitZIG_two_step(Y_roi, age)
% FITZIG_TWO_STEP - Two-step estimation of Zero-Inflated Gaussian (ZIG) model
%
% INPUTS:
%   Y_roi : [Nv x N] matrix of voxel-level data (voxels × subjects)
%   age   : [N x 1] vector of subject ages
%
% OUTPUT:
%   result : struct with fields:
%     - beta     : [Nv x 3] regression coefficients for Gaussian part
%     - gamma    : [Nv x 2] regression coefficients for zero inflation
%     - sigma    : [Nv x 1] residual standard deviation (Gaussian)
%     - loglik   : [Nv x 1] approximate log-likelihood
%     - iter     : [Nv x 1] number of iterations (always 1 for 2-step)
%     - p_beta   : [Nv x 3] p-values for Gaussian coefficients
%     - p_gamma  : [Nv x 2] p-values for logistic coefficients
%     - mu_est   : [Nv x N] conditional mean predictions (only where Y > 0)
%     - pi_est   : [Nv x N] predicted zero-inflation probabilities

    [Nv, N] = size(Y_roi);
    age = age(:);  % Ensure column vector

    % Standardize age
    age_std = (age - mean(age)) / std(age);

    % Design matrices
    X_mu = [ones(N,1), age_std, age_std.^2];  % [N x 3]
    X_pi = [ones(N,1), age_std];              % [N x 2]

    % Preallocate outputs
    beta     = zeros(Nv, 3);
    gamma    = zeros(Nv, 2);
    sigma    = zeros(Nv, 1);
    loglik   = zeros(Nv, 1);
    iter     = ones(Nv, 1);
    p_beta   = zeros(Nv, 3);
    p_gamma  = zeros(Nv, 2);
    mu_est   = zeros(Nv, N);
    pi_est   = zeros(Nv, N);

    % Options for logistic regression
    opts = statset('fitglm');
    opts.MaxIter = 1000;

    for v = 1:Nv
        y = Y_roi(v,:)';
        is_pos = y > 0;

        %% === STEP 1: Logistic Regression ===
        if all(is_pos) || all(~is_pos)
            gamma(v,:)   = 0;
            p_gamma(v,:) = 1;
            pi_est(v,:)  = mean(~is_pos) * ones(1,N);
        else
            try
                glm = fitglm(X_pi, is_pos, ...
                    'Distribution', 'binomial', ...
                    'Link', 'logit', ...
                    'Intercept', false, ...
                    'Options', opts);
                gamma(v,:)   = glm.Coefficients.Estimate';
                p_gamma(v,:) = glm.Coefficients.pValue';
                pi_hat       = glm.Fitted.Probability';  % Pr(Y > 0)
                pi_est(v,:)  = 1 - pi_hat;               % p = Pr(Y = 0)
            catch
                gamma(v,:)   = 0;
                p_gamma(v,:) = 1;
                pi_est(v,:)  = mean(~is_pos) * ones(1,N);
            end
        end

        %% === STEP 2: Gaussian Regression on Positives ===
        if sum(is_pos) > 10  % Conservative threshold
            try
                X_mu_pos = X_mu(is_pos,:);
                y_pos    = y(is_pos);

                % Check for rank deficiency
                if rank(X_mu_pos) < size(X_mu_pos,2)
                    warning('Rank deficiency at voxel %d; skipping linear regression.', v);
                    continue;
                end

                lm = fitlm(X_mu_pos, y_pos, 'Intercept', false);
                beta(v,:)   = lm.Coefficients.Estimate';
                p_beta(v,:) = lm.Coefficients.pValue';
                sigma(v)    = sqrt(lm.MSE);

                % Predict conditional means
                mu_pos = X_mu_pos * beta(v,:)';
                mu_est(v, is_pos) = mu_pos;

                % Approximate log-likelihood
                pi_v = pi_est(v,:);
                ll_zeros = sum(log(pi_v(~is_pos) + eps));
                ll_pos   = sum(log(1 - pi_v(is_pos) + eps)) + ...
                           sum(log(normpdf(y_pos, mu_pos, sigma(v)) + eps));
                loglik(v) = ll_zeros + ll_pos;
            catch
                beta(v,:)   = 0;
                p_beta(v,:) = 1;
                sigma(v)    = 0;
            end
        end
    end

    % Assemble result
    result.beta    = beta;
    result.gamma   = gamma;
    result.sigma   = sigma;
    result.loglik  = loglik;
    result.iter    = iter;
    result.p_beta  = p_beta;
    result.p_gamma = p_gamma;
    result.mu_est  = mu_est;
    result.pi_est  = pi_est;
end

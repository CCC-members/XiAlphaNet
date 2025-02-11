function k_min = findMinimumK(parameters, RE_tol, Nxrand)
    % Parameters:
    % parameters: structure containing all required settings
    % RE_tol: Relative Error percentage
    % Nxrand: number of random x trials

    Nv = parameters.Dimensions.Nv;
    xdim = 7 * Nv + 1;

    % Initialize variables
    max_k = length(parameters.Data.freq); % Upper bound for k to prevent infinite loops
    relative_error_threshold = RE_tol / 100; % Error < relative error threshold
    k_min = max_k; % Start with the largest k and decrease if possible

    % Preallocate array to store mean relative errors for each k
    mean_relative_errors = zeros(1, max_k);

    % Generate a random x sample for the entire computation
    x = generateRandomSample(parameters, 10); % Generate random x
    x,Ne,T,sw,sp,nsf_band,Sw
    % Evaluate F fully without stochastic sampling
    parameters.Stochastic.stoch = 1;
    parameters.Stochastic.Nsfreq = 47; % Set current k
    parameters = sample_frequencies(parameters);
    F_full = evaluateF(x, parameters);
    [dF_full, ~, ~] = evaluatedF(x, parameters);

    % Outer loop over k
    tic
    for k = 1:max_k
        relative_errors = zeros(1, Nxrand); % Track relative errors for current k

        % Parallelized inner loop for Nxrand random trials
        for trial = 1:Nxrand

            % Adjust parameters for stochastic evaluation
            param_trial = parameters; % Create a local copy for parallel safety
            param_trial.Stochastic.stoch = 1; % Enable stochastic evaluation
            param_trial.Stochastic.Nsfreq = k; % Set current k
            param_trial = sample_frequencies(param_trial);

            % Stochastic evaluation of F
            F_stoch = evaluateF(x, param_trial);
            [dF_stoch, ~, ~] = evaluatedF(x, param_trial);

            % Calculate relative absolute error
            if norm(F_full) ~= 0 && norm(dF_full) ~= 0
                relative_error = abs(F_full - F_stoch) / abs(F_full);
                df_relative_error = norm(dF_full - dF_stoch, 'fro') / norm(dF_full, 'fro');
            else
                % Handle the case where F_full or dF_full is zero to avoid division by zero
                relative_error = abs(F_full - F_stoch) / 1e-10;
                df_relative_error = norm(dF_full - dF_stoch, 'fro') / 1e-10;
            end

            % Average the relative errors
            relative_error = (relative_error + df_relative_error) / 2;

            % Store the relative error for this trial
            relative_errors(trial) = relative_error;
        end

        % Calculate and store the mean relative error for this k
        mean_relative_errors(k) = mean(relative_errors);

        % Check if the mean relative error for this k is within the threshold
        if mean_relative_errors(k) < relative_error_threshold
            k_min = k; % Update k_min if the new k is more efficient
            break; % Stop searching as we have found the minimum k satisfying the condition
        end
    end

    % Output the minimum k or an error if no solution is found within bounds
    if k_min == max_k
        error('Failed to achieve the desired accuracy with k up to %d\n', max_k);
    endlambda2
toc
    % Optional: Plot the mean relative errors for each k (can be uncommented)
    % figure;
    % plot(1:max_k, movmean(mean_relative_errors * 100, 1), '-o');
    % xlabel('k');
    % ylabel('Mean Relative Error (%)');
    % title('Mean Relative Error across Simulations for each k');
    % grid on;
end
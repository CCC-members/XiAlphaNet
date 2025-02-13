function k_min = findMinimumK(freq, T, Sw, RE_tol, Nxrand, index_parallel)
    % Parameters:
    % RE_tol: Relative Error percentage
    % Nxrand: number of random x trials
    % index_parallel: Whether to use parallel processing (1: parallel, 0: serial)

    [Ne, Nv, ~] = size(T);
    xdim = 7 * Nv + 1;

    % Initialize variables
    max_k = length(freq); % Upper bound for k to prevent infinite loops
    relative_error_threshold = RE_tol / 100; % Error < relative error threshold
    k_min = max_k; % Start with the largest k and decrease if possible
    % Preallocate array to store mean relative errors for each k
    mean_relative_errors = zeros(1, max_k);

    % Generate a random x sample for the entire computation
    x = generateRandomSample(Ne, Nv, Sw, T, freq, 10); % Generate random x
    % Evaluate F fully without stochastic sampling
    index_stoch = 0;
    Nsfreq = 47; % Set current k
    [nsf_band, sw, sp] = sample_frequencies(freq, index_stoch, Nsfreq);
    F_full = evaluateF(x, Ne, T, sw, sp, nsf_band, Sw);
    [dF_full, ~, ~] = evaluatedF(x, Ne, Nv, T, sw, sp, nsf_band, Sw);

    % Outer loop over k
    tic
    for k = 1:max_k
        relative_errors = zeros(1, Nxrand); % Track relative errors for current k

        if index_parallel == 1
            % Parallelized inner loop for Nxrand random trials
            parfor trial = 1:Nxrand
                % Adjust parameters for stochastic evaluation
                index_stoch = 1;
                Nsfreq_trial = k; % Set current k
                [nsf_band_trial, sw_trial, sp_trial] = sample_frequencies(freq, index_stoch, Nsfreq_trial);

                % Stochastic evaluation of F
                F_stoch = evaluateF(x, Ne, T, sw_trial, sp_trial, nsf_band_trial, Sw);
                [dF_stoch, ~, ~] = evaluatedF(x, Ne, Nv, T, sw_trial, sp_trial, nsf_band_trial, Sw);

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
        else
            % If parallelization is not used, use a regular for loop
            for trial = 1:Nxrand
                % Adjust parameters for stochastic evaluation
                index_stoch = 1;
                Nsfreq_trial = k; % Set current k
                [nsf_band_trial, sw_trial, sp_trial] = sample_frequencies(freq, index_stoch, Nsfreq_trial);

                % Stochastic evaluation of F
                F_stoch = evaluateF(x, Ne, T, sw_trial, sp_trial, nsf_band_trial, Sw);
                [dF_stoch, ~, ~] = evaluatedF(x, Ne, Nv, T, sw_trial, sp_trial, nsf_band_trial, Sw);

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
        end

        % Calculate and store the mean relative error for this k
        mean_relative_errors(k) = mean(relative_errors);

        % Check if the mean relative error for this k is within the threshold
        if mean_relative_errors(k) < relative_error_threshold
            k_min = k; % Update k_min if the new k is more efficient
            break; % Stop searching as we have found the minimum k satisfying the condition
        end
    end

    % % Output the minimum k or an error if no solution is found within bounds
    % if k_min == max_k
    %     error('Failed to achieve the desired accuracy with k up to %d\n', max_k);
    % end
    toc
end

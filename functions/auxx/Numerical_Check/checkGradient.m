function checkGradient(parameters, K, T)
    % K: number of samples
    % T: number of target components to evaluate stochastically


    % Assumption of dimensions, Nr and Ne should be defined or passed as needed
    Nr = parameters.Dimensions.Nr;
    Ne = parameters.Dimensions.Ne;
    parameters = sample_frequencies(parameters);
    
    if nargin<3
        K = 200;
        T = 1;
    end


    % Initialize accumulators for the errors
    total_relative_error = 0;

    for k = 1:K
        % Generate random sample for x
        e = rand(Nr, 3); % Example matrix e
        a = rand(Nr, 4); % Example matrix a
        sigma2 = rand(); % Example sigma^2
        x = [e(:);a(:);sigma2]; % Concatenate into vector

        % Compute the analytical gradient
        [dF_analytical, ~,~] = evaluatedF(x, parameters);

        % Compute the numerical gradient for a stochastic subset of components
        epsilon = 1e-4; % Perturbation size
        dF_numerical = zeros(length(x), 1);
        indices = randperm(length(x), T); % Randomly select T indices to evaluate

        for i = indices
            x_plus = x;
            x_minus = x;
            x_plus(i) = x_plus(i) + epsilon;
            x_minus(i) = x_minus(i) - epsilon;
            tic;
            sF_plus = evaluateF(x_plus, parameters);
            t = toc;
            sF_minus = evaluateF(x_minus, parameters);
            dF_numerical(i) = (sF_plus - sF_minus) / (2 * epsilon);
            fprintf('--> Partial derivative Position %d, Simulation %d, Time %.4f seconds\n', i, k, t);
        end

        % Compute the relative absolute error for sampled indices
        relative_error = compute_relative_error(dF_analytical, dF_numerical, indices);

        % Accumulate the total relative error
        total_relative_error = total_relative_error + relative_error;
    end

    % Compute the mean relative absolute error
    mean_relative_error = total_relative_error / K;

% Display the results
fprintf('Mean relative absolute error between analytical and numerical gradients: %e\n', mean_relative_error*100);

% Check if mean relative error exceeds 10%
if mean_relative_error > 0.1
    warning('Mean relative error exceeds 1%%: %e\n', mean_relative_error*100);
end

figure(1)
hold on
plot(dF_analytical(indices), 'o', 'Color', 'b')
plot(dF_numerical(indices), '*', 'Color', 'r')

function relative_error = compute_relative_error(dF_analytical, dF_numerical, indices)
    % Compute the Frobenius norms
    norm_num = norm(dF_numerical(indices), 'fro');
    norm_difference = norm(dF_analytical(indices) - dF_numerical(indices), 'fro');
    
    % Avoid division by zero
    if norm_num == 0
        norm_num = 1e-20;
    end
    
    % Compute the relative error
    relative_error = norm_difference / norm_num;
end


end



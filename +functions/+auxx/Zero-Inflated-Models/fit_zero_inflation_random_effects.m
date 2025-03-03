function [model_output, temp_model] = fit_zero_inflation_random_effects(x, y, group, x_eval)
    % fit_zero_inflation_random_effects - Fit a zero-inflated model with random effects using R and obtain predictions and SEs for x_eval.
    %
    % Inputs:
    %   x        - Predictor variable, a numeric vector.
    %   y        - Response variable (count data), a numeric vector.
    %   group    - Grouping variable for random effects (e.g., subject IDs).
    %   x_eval   - New points at which to evaluate the model (numeric vector).
    %
    % Outputs:
    %   model_output - A structure containing:
    %       1. model coefficients (conditional and zero-inflation).
    %       2. predictions for both conditional and zero-inflation models at x_eval.
    %       3. standard errors for both models at x_eval.

    % Write the data to a temporary CSV file for R
    if nargin < 4 || isempty(x_eval)
        error('x_eval is required for predictions.');
    end
    
    if nargin < 3 || isempty(group) % If no group is provided, do not use random effects
        group = []; 
        T = table(x(:), y(:), 'VariableNames', {'x', 'y'}); % Without group column
    else
        T = table(x(:), y(:), group(:), 'VariableNames', {'x', 'y', 'group'}); % With group column
    end
    writetable(T, 'temp_data.csv');
    
    % Pass x_eval to a CSV file for R evaluation
    x_eval_table = table(x_eval(:), 'VariableNames', {'x'});
    writetable(x_eval_table, 'x_eval.csv');

    % Run the R script from MATLAB
    [status, cmdout] = system('"C:\Program Files\R\R-4.4.2\bin\Rscript" fit_zero_inflation_model2.R');
    
    % Check for successful execution
    if status ~= 0
        error('Error running R script: %s', cmdout);
    end
    
    % Load the model coefficients and predictions from CSV files
    temp_model = readtable('model_coefficients.csv');
    predictions = readtable('predictions_with_se.csv');
    
    % Clean up temporary files
    delete('temp_data.csv');
    delete('x_eval.csv');
    delete('model_coefficients.csv');
    delete('predictions_with_se.csv');
    
    % Convert table to structure for easier access in MATLAB
    temp_model = table2struct(temp_model);
    predictions = table2struct(predictions);
    
    % Extract coefficients (conditional and zero-inflation)
    cond_coefs = temp_model(1:3);
    zi_coefs = temp_model(4:5); %6
    
    % Convert coefficients (estimate values) to numerical arrays for further calculation
    beta_conditional = cell2mat({cond_coefs.estimate});
    beta_zero_inflation = cell2mat({zi_coefs.estimate});
    
    % Initialize the output structure for model details
    model_output = struct();
    
    % Store coefficients in the model_output structure
    model_output.coef.cond = beta_conditional;
    model_output.coef.zi = beta_zero_inflation;
    
    % Store predictions and standard errors
    model_output.predictions = struct();
    model_output.predictions.cond = [predictions.cond_pred];
    model_output.predictions.cond_se = [predictions.cond_se];
    model_output.predictions.zi = [predictions.zi_pred];
    model_output.predictions.zi_se = [predictions.zi_se];
end

% %% example
% % Generate synthetic data with n = 50
% n = 50;  % Number of data points
% x = linspace(0, 10, n);  % Predictor variable (e.g., time or feature)
% group = randi([1, 5], n, 1);  % Random grouping variable (5 groups)
% beta_conditional = [-0.1, 0.3, -0.02];  % Coefficients for conditional model
% beta_zero_inflation = [2, -1.5, 0.1];  % Coefficients for zero-inflation model
% 
% % Generate conditional model (log(mu) = beta0 + beta1*x + beta2*x^2)
% mu = exp(beta_conditional(1) + beta_conditional(2)*x + beta_conditional(3)*x.^2);
% 
% % Generate the zero-inflation model (logit(pi) = beta0 + beta1*x + beta2*x^2)
% prob_zero = 1 ./ (1 + exp(-(beta_zero_inflation(1) + beta_zero_inflation(2)*x + beta_zero_inflation(3)*x.^2)));
% 
% % Generate response variable with random zeros based on zero-inflation probability
% y = zeros(1, n);
% for i = 1:n
%     if rand < prob_zero(i)  % With probability prob_zero, set y to 0
%         y(i) = 0;
%     else
%         y(i) = poissrnd(mu(i));  % Otherwise, generate from Poisson distribution
%     end
% end
% 
% % Create a table for data
% data = table(x', y', group, 'VariableNames', {'x', 'y', 'group'});
% 
% model_output = fit_zero_inflation_random_effects(data.x, data.y, data.group, x_eval);
% cond_pred = model_output.predictions.cond;  % Conditional model predictions
% cond_se = model_output.predictions.cond_se;  % Conditional model standard errors
% zi_pred = model_output.predictions.zi;  % Zero-inflation model predictions
% zi_se = model_output.predictions.zi_se;  % Zero-inflation model standard errors
% 
% % Calculate the upper and lower bounds for the conditional model (prediction +/- error)
% cond_upper = cond_pred + cond_se;  % Upper bound for conditional model (prediction + error)
% cond_lower = cond_pred - cond_se;  % Lower bound for conditional model (prediction - error)
% 
% % Calculate the upper and lower bounds for the zero-inflation model (prediction +/- error)
% zi_upper = zi_pred + zi_se;  % Upper bound for zero-inflation model (prediction + error)
% zi_lower = zi_pred - zi_se;  % Lower bound for zero-inflation model (prediction - error)
% 
% % Create the figure
% figure;
% 
% % Plot the conditional model with error bands
% subplot(2, 1, 1);  % Create a subplot for conditional model
% hold on;
% plot(x_eval, cond_pred, 'b', 'LineWidth', 2);  % Plot the conditional predictions (central curve)
% plot(x_eval, cond_upper, 'r--', 'LineWidth', 1);  % Plot the upper error bound
% plot(x_eval, cond_lower, 'r--', 'LineWidth', 1);  % Plot the lower error bound
% xlabel('x_eval');
% ylabel('Conditional Model Prediction');
% title('Conditional Model with Error Estimator');
% legend('Conditional Prediction', 'Upper Error Bound', 'Lower Error Bound', 'Location', 'Best');
% hold off;
% 
% % Plot the zero-inflation model with error bands
% subplot(2, 1, 2);  % Create a subplot for zero-inflation model
% hold on;
% plot(x_eval, zi_pred, 'g', 'LineWidth', 2);  % Plot the zero-inflation predictions (central curve)
% plot(x_eval, zi_upper, 'm--', 'LineWidth', 1);  % Plot the upper error bound
% plot(x_eval, zi_lower, 'm--', 'LineWidth', 1);  % Plot the lower error bound
% xlabel('x_eval');
% ylabel('Zero-Inflation Model Prediction');
% title('Zero-Inflation Model with Error Estimator');
% legend('Zero-Inflation Prediction', 'Upper Error Bound', 'Lower Error Bound', 'Location', 'Best');
% hold off;

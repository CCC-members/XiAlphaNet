function model_output = fit_zero_inflated_model(x, y, x_eval)
    % fit_zero_inflated_model - Fit and evaluate a zero-inflated Poisson model.
    %
    % This function fits a zero-inflated Poisson model to count data `y` 
    % with predictor variable `x`, using Python's `statsmodels` library. 
    % Zero-inflated models are useful when count data contains an excess 
    % of zeros that cannot be explained by a standard Poisson or 
    % negative binomial distribution alone.
    %
    % A zero-inflated Poisson model consists of two parts:
    %   1. A count model (e.g., Poisson) that explains the distribution 
    %      of non-zero counts.
    %   2. A zero-inflation model, often with a logistic regression, 
    %      that explains the probability of generating an excess zero.
    %
    % This function integrates MATLAB with Python by calling a custom 
    % Python function `fit_zero_inflated_model` (in 'zero_inflation_model.py') 
    % from MATLAB using Python's `statsmodels` library.
    %
    % Inputs:
    %   x - Predictor variable, a numeric vector.
    %   y - Response variable (count data), a numeric vector.
    %   x_eval - (Optional) A vector of evaluation points for the predictor
    %            variable, at which the model will be evaluated.
    %
    % Outputs:
    %   model_output - A structure containing model coefficients, p-values,
    %                  fit statistics (AIC, BIC), and, if `x_eval` is provided,
    %                  the predicted counts and zero-inflation probabilities at `x_eval`.
    %
    % Requirements:
    %   - Python 3.9 or compatible version with `statsmodels` installed.
    %   - `zero_inflation_model.py` script in the MATLAB path.

    % Convert MATLAB data to Python lists
    py_y = py.list(y);  % Convert y to a Python list
    py_x = py.list(x);  % Convert x to a Python list

    % Import the Python module (assuming 'zero_inflation_model.py' is in the path)
    pyModule = py.importlib.import_module('zero_inflation_model');

    % Call the Python function to fit the zero-inflated model
    result = pyModule.fit_zero_inflated_model(py_y, py_x);

    % Convert the Python dictionary result to a MATLAB structure
    output = struct(result);

    % Create a MATLAB structure to hold the outputs
    model_output = struct();

    % Convert and store each part in the output structure
    model_output.params = struct(output.params);   % Model coefficients
    model_output.pvalues = struct(output.pvalues); % P-values for each coefficient
    model_output.aic = double(output.aic);         % AIC value
    model_output.bic = double(output.bic);         % BIC value

    % If x_eval is provided, evaluate the model at the specified predictor values
    if nargin > 2
        % Extract model parameters
        intercept = model_output.params.const;          % Intercept for count model
        beta_x = model_output.params.x;                 % Coefficient for x in count model
        beta_x2 = model_output.params.x2;               % Coefficient for x^2 in count model
        zi_intercept = model_output.params.inflate_const; % Intercept for zero-inflation model
        zi_beta_x = model_output.params.inflate_x;      % Coefficient for x in zero-inflation model
        zi_beta_x2 = model_output.params.inflate_x2;    % Coefficient for x^2 in zero-inflation model

        % Calculate predicted counts using the count model
        predicted_counts = exp(intercept + beta_x * x_eval + beta_x2 * (x_eval .^ 2));

        % Calculate zero-inflation probability using the logit link
        zi_prob = 1 ./ (1 + exp(-(zi_intercept + zi_beta_x * x_eval + zi_beta_x2 * (x_eval .^ 2))));

        % Store evaluations in the output structure
        model_output.x_eval = x_eval;                  % Evaluation points for the predictor
        model_output.predicted_counts = predicted_counts; % Predicted count values at x_eval
        model_output.zi_prob = zi_prob;                % Zero-inflation probabilities at x_eval
    end
end

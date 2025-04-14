function [lambda] = lambda_regspace(freq, T, Sw, index_stoch, Nsfreq, x0)
    import functions.auxx.StochasticEval.*
    import functions.auxx.GenerateSourceSample.*
    import functions.FunctionGrandientProx.*
    import functions.auxx.RegSpace.*
    import functions.auxx.ModelVectorization.*

    [Ne, Nv, ~] = size(T);
    [nsf_band, sw, sp] = sample_frequencies(freq, index_stoch, Nsfreq);

    % Generate a dense log-spaced range of epsilon values
    epsilon_values = logspace(-2, -40, 20);  % 20 values from 1e-2 to 1e-20

    % Initialize lambda estimates to +Inf (for minimum selection)
    lambda_best = [Inf, Inf, Inf];

    for eps = epsilon_values
        dim = size(x0);
        x_init = x0;

        % Evaluate gradient at this scale
        dF = evaluatedF(x_init, Ne, Nv, T, sw, sp, nsf_band, Sw);
        [dFe, dFa, dFs2] = x2v(dF);

        % Compute lambdas for each block
        lambda_s2 = compute_lambda(dFs2)^2;  % Ridge block
        lambda_e  = compute_lambda(dFe);     % Group Lasso
        lambda_a  = compute_lambda(dFa);     % Group Lasso

        % Retain the minimum (most conservative) across all epsilons
        lambda_best(1) = min(lambda_best(1), lambda_s2);
        lambda_best(2) = min(lambda_best(2), lambda_e);
        lambda_best(3) = min(lambda_best(3), lambda_a);
    end

    % Apply safety margin to avoid null solution
    lambda = lambda_best * 0.9;
end


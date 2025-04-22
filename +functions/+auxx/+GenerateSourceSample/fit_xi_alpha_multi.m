function [best_params, best_fit, all_fits] = fit_xi_alpha_multi(S, omega, N_init, doplot)
    % Fit Xi-Alpha model to a single spectrum using multiple initializations
    % Inputs:
    %   S       - spectrum vector (Nw x 1), real-valued
    %   omega   - frequency vector (Nw x 1)
    %   N_init  - number of initializations to try
    %   doplot  - (optional) if true, plot in dB scale
    %
    % Outputs:
    %   best_params - best parameter estimate [e1, e2, e3, a1, a2, a3, a4]
    %   best_fit    - fitted spectrum using best_params
    %   all_fits    - cell array of all fitted spectra

    if nargin < 4
        doplot = false;
    end

    % Ensure input is real-valued and column vectors
    S = real(S(:));  
    omega = omega(:);

    % === Xi and Alpha model definitions ===
    xi_omega = @(e, omega) e(1) ./ (1 + e(2) * omega.^2).^e(3);
    alpha_omega = @(a, omega) a(1) ./ (1 + a(2) * (omega - a(4)).^2).^a(3);
    model_fun = @(params, omega) ...
        xi_omega(params(1:3), omega) + alpha_omega(params(4:7), omega);

    % === Bounds ===
    lb = [0, 0, 0,     0, 0, 0, 7];     % [e1, e2, e3, a1, a2, a3, a4]
    ub = [Inf, 10, 10,  Inf, 10, 100, 13];

    % === Optimizer options ===
    options = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxIterations', 1000, 'FunctionTolerance', 1e-8);

    % === Preallocation ===
    all_params = zeros(N_init, 7);
    all_resnorm = zeros(N_init, 1);
    all_fits = cell(N_init, 1);

    % === Fix e1 from the lowest-frequency value ===
    e1_fixed = S(1);

    % === Run multiple initializations ===
    for k = 1:N_init
        % Generate random initial parameters
        e2_0 = rand() * 1.0;
        e3_0 = 1 + rand() * 1.5;
        a1_0 = max(S) * (0.5 + rand());
        a2_0 = rand() * 0.5;
        a3_0 = 2 + rand() * 3;
        a4_0 = 10 + rand() * 2;

        init_params = [e1_fixed, e2_0, e3_0, a1_0, a2_0, a3_0, a4_0];

        % Fit
        [params_k, resnorm_k] = lsqcurvefit( ...
            @(params, omega) model_fun(params, omega), ...
            init_params, omega, S, lb, ub, options);

        all_params(k, :) = params_k;
        all_resnorm(k) = resnorm_k;
        all_fits{k} = model_fun(params_k, omega);
    end

    % === Select best result ===
    [~, best_idx] = min(all_resnorm);
    best_params = all_params(best_idx, :);
    best_fit = all_fits{best_idx};

    % === Optional Plotting ===
    if doplot
        figure;
        plot(omega, 10*log10(S), 'k.-', 'DisplayName', 'Observed');
        hold on;
        plot(omega, 10*log10(best_fit), 'r-', 'LineWidth', 2, 'DisplayName', 'Xi-Alpha Fit');
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        title('Xi-Alpha Fit (Log Scale)');
        legend('show');
        grid on;
    end
end

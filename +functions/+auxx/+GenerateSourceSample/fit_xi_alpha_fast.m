function [params, fit] = fit_xi_alpha_fast_log(S, omega, doplot)
    % Fast Xi-Alpha fitting in log-space
    % Inputs:
    %   S       - spectrum (Nw x 1), positive and real
    %   omega   - frequency vector (Nw x 1)
    %   doplot  - optional flag to enable plotting
    %
    % Outputs:
    %   params  - [e1, e2, e3, a1, a2, a3, a4]
    %   fit     - fitted spectrum (linear scale)

    if nargin < 3
        doplot = false;
    end

    % Ensure column vectors and positivity
    S = max(real(S(:)), eps);  % Avoid log(0)
    omega = omega(:);

    % === Fixed parameters ===
    e1 = S(1);  % Low-frequency offset

    % Restrict to alpha band (713 Hz)
    alpha_idx = (omega >= 7) & (omega <= 13);
    omega_alpha = omega(alpha_idx);
    S_alpha = S(alpha_idx);
    [a1, idx_max] = max(S_alpha);
    a4 = omega_alpha(idx_max);

    % === Model in log-space ===
    model_log = @(p, omega) log( ...
        e1 ./ (1 + p(1)*omega.^2).^p(2) + ...
        a1 ./ (1 + p(3)*(omega - a4).^2).^p(4) );

    % === Optimization ===
    p0 = [0.1, 1.5, 0.1, 3];      % [e2, e3, a2, a3]
    lb = [0,   0.5, 0,   1];
    ub = [10, 10,  2,  10];

    options = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxIterations', 500, 'FunctionTolerance', 1e-6);

    logS = log(S);
    p_opt = lsqcurvefit(@(p, omega) model_log(p, omega), ...
        p0, omega, logS, lb, ub, options);

    % === Final model evaluation (in linear scale) ===
    fit = e1 ./ (1 + p_opt(1)*omega.^2).^p_opt(2) + ...
          a1 ./ (1 + p_opt(3)*(omega - a4).^2).^p_opt(4);

    % Output all parameters
    params = [e1, p_opt(1:2), a1, p_opt(3:4), a4];

    % === Plot ===
    if doplot
        figure;
        plot(omega, 10*log10(S), 'k.-', 'DisplayName', 'Observed');
        hold on;
        plot(omega, 10*log10(fit), 'r-', 'LineWidth', 2, 'DisplayName', 'Xi-Alpha Log Fit');
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        title('Xi-Alpha Model Fit in Log Scale');
        legend('show'); grid on;
    end
end

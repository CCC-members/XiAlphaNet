function fooof_results = fooof_matlab(freqs, power_spectrum, freq_range, n_peaks)
% FOOOF_MATLAB: MATLAB-native implementation of FOOOF (Fitting Oscillations and One-Over-F)
%
% Inputs:
%   freqs            - Frequency vector [Hz]
%   power_spectrum   - Power spectrum [same length as freqs]
%   freq_range       - [low_freq, high_freq] analysis range [Hz]
%   n_peaks          - Maximum number of peaks to fit (default: 3)
%
% Outputs:
%   fooof_results    - Struct containing:
%       .aperiodic_params - [offset, slope]
%       .peak_params      - [amplitude, center_freq, width] for each peak
%       .model_fit        - Reconstructed model fit in log10 power
%       .r_squared        - Goodness-of-fit metric

    if nargin < 4
        n_peaks = 3;  % Default to 3 peaks
    end

    % Subset the frequency range
    idx_range = (freqs >= freq_range(1)) & (freqs <= freq_range(2));
    freqs_fit = freqs(idx_range);
    power_fit = power_spectrum(idx_range);

    % Log-transform the data (log10 power)
    log_freqs = log10(freqs_fit(:));  % ensure column vector
    log_power = log10(abs(power_fit(:)));  % ensure column vector

    % --- Fit the aperiodic (1/f) component ---
    X = [ones(length(log_freqs),1), -log_freqs];  % offset and negative slope
    b = real(X \ log_power);  % Least squares linear fit, ensure real

    offset = real(b(1));
    slope  = real(b(2));

    % Remove the aperiodic component
    aperiodic_fit = X * b;
    residual = log_power - aperiodic_fit;

    % --- Fit peaks (Gaussian) in the residual ---
    peak_params = [];
    residual_fit = zeros(size(residual));

    for peak_idx = 1:n_peaks
        % Find peak in residual
        [amp, idx_peak] = max(residual);
        if amp < 0.05  % Stop if no significant peak
            break;
        end
        center_freq = freqs_fit(idx_peak);

        % Estimate width by finding half-maximum points
        half_max = amp / 2;

        % LEFT side
        left_idx = find(residual(1:idx_peak) < half_max, 1, 'last');
        if isempty(left_idx)
            left_idx = 1;
        end
        left_freq = freqs_fit(left_idx);

        % RIGHT side
        right_idx_relative = find(residual(idx_peak:end) < half_max, 1, 'first');
        if isempty(right_idx_relative)
            right_idx = length(freqs_fit);
        else
            right_idx = idx_peak - 1 + right_idx_relative;
            right_idx = min(right_idx, length(freqs_fit));
        end
        right_freq = freqs_fit(right_idx);

        % Width estimate (standard deviation)
        est_width = max(0.5, (right_freq - left_freq)/2.355);  % approximate sigma

        % Define the Gaussian in log-power space
        gauss_func = @(p, x) real(p(1) * exp(-0.5 * ((x(:) - p(2)) ./ p(3)).^2));

        % Initial guess and bounds
        p0 = [amp, center_freq, est_width];
        lb = [0, center_freq - 2, 0.1];
        ub = [Inf, center_freq + 2, 10];
        options = optimset('Display', 'off');

        % Fit the Gaussian
        [p_fit, ~] = lsqcurvefit(gauss_func, p0, freqs_fit(:), residual, lb, ub, options);
        p_fit = real(p_fit);

        % Store the peak
        peak_params = [peak_params; p_fit];

        % Subtract the fitted peak from the residual
        residual = residual - gauss_func(p_fit, freqs_fit(:));
        residual_fit = residual_fit + gauss_func(p_fit, freqs_fit(:));
    end

    % --- Reconstruct the full model in log10 power ---
    model_fit = aperiodic_fit + residual_fit;

    % --- Goodness-of-fit ---
    r_squared = 1 - sum((log_power - model_fit).^2) / sum((log_power - mean(log_power)).^2);

    % --- Pack results ---
    fooof_results = struct();
    fooof_results.aperiodic_params = [offset, slope];
    fooof_results.peak_params = peak_params; % each row: [amplitude, center_freq, width]
    fooof_results.model_fit = model_fit;
    fooof_results.r_squared = real(r_squared);

    % % --- Plot results ---
    % figure;
    % plot(freqs_fit, log_power, 'k', 'LineWidth', 1.5); hold on;
    % plot(freqs_fit, model_fit, 'r', 'LineWidth', 2);
    % xlabel('Frequency (Hz)');
    % ylabel('Log10 Power');
    % legend('Original Spectrum', 'FOOOF Fit');
    % title('FOOOF MATLAB Implementation (Log10 Power)');
    % grid on;

end

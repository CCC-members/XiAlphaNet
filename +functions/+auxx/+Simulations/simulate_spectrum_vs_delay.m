function simulate_spectrum_vs_delay(C, varargin)
% SIMULATE_SPECTRUM_VS_DELAY Visualizes how delays affect spectral power
% 
% Usage:
%   simulate_spectrum_vs_delay(C)
%   simulate_spectrum_vs_delay(C, 'ref_voxel', 30, 'freqs', 1:0.5:40, ...)
%
% Inputs:
%   C           - [Nv x Nv] connectivity matrix (must be square, zero-diagonal)
%
% Optional Name-Value Pairs:
%   'ref_voxel' - Index of reference voxel (default = center)
%   'freqs'     - Frequency range (Hz) (default = linspace(1, 40, 200))
%   'delays'    - Delay range (ms) (default = linspace(6, 30, 100))
%
% Output: 4-view plot of log spectral power vs. frequency and delay

% Parse inputs
Nv = size(C,1);
p = inputParser;
addParameter(p, 'ref_voxel', round(Nv/2));
addParameter(p, 'freqs', linspace(1, 40, 200));
addParameter(p, 'delays', linspace(6, 30, 100));
parse(p, varargin{:});

ref_voxel = p.Results.ref_voxel;
freqs = p.Results.freqs(:);
delays = p.Results.delays(:);
omega = 2 * pi * freqs;
Nf = numel(freqs);
Nd = numel(delays);

% Alpha frequency shift params
F_alpha_0 = 10;      % Hz at minimum delay
delta_F = 4;         % total shift (Hz)

% Fixed xi power (sharper decay)
xi_power = (1 ./ (1 + (omega / (2*pi*7)).^2)).^3.5;

% Initialize
log_power_map = zeros(Nf, Nd);
I = eye(Nv);

for di = 1:Nd
    d = delays(di);  % current delay (ms)

    % Delay-dependent alpha peak frequency
    shift_ratio = (d - min(delays)) / (max(delays) - min(delays));
    F_alpha_d = F_alpha_0 - delta_F * shift_ratio;

    alpha_power = (1 ./ (1 + ((omega - 2*pi*F_alpha_d) / (2*pi*5)).^2)).^2;
    total_power = repmat(xi_power + alpha_power, 1, Nv);

    for fi = 1:Nf
        w = omega(fi);
        A = zeros(Nv);
        A(ref_voxel, :) = C(ref_voxel, :) * exp(-1i * w * d);
        A(ref_voxel, ref_voxel) = 0;

        try
            G = inv(I - A);
        catch
            G = pinv(I - A);  % fallback
        end

        g_row = G(ref_voxel, :);
        S = sum(abs(g_row).^2 .* total_power(fi, :));
        log_power_map(fi, di) = 10 * log10(S + 1e-10);
    end
end

% Plot
[DelayGrid, FreqGrid] = meshgrid(delays, freqs);

figure('Position', [100 100 1400 800]);

subplot(2,2,1);
surf(DelayGrid, FreqGrid, log_power_map, 'EdgeColor', 'none');
view(0, 90); xlabel('Delay (ms)'); ylabel('Frequency (Hz)');
zlabel('Log Power (dB)'); title('Top View'); colormap parula;

subplot(2,2,2);
surf(DelayGrid, FreqGrid, log_power_map, 'EdgeColor', 'none');
view(45, 30); xlabel('Delay (ms)'); ylabel('Frequency (Hz)');
zlabel('Log Power (dB)'); title('Side View'); colormap parula;

subplot(2,2,3);
surf(DelayGrid, FreqGrid, log_power_map, 'EdgeColor', 'none');
view(120, 45); xlabel('Delay (ms)'); ylabel('Frequency (Hz)');
zlabel('Log Power (dB)'); title('Oblique View'); colormap parula;

subplot(2,2,4);
imagesc(delays, freqs, log_power_map); axis xy;
xlabel('Delay (ms)'); ylabel('Frequency (Hz)');
title('Heatmap View'); colormap parula;
c = colorbar; ylabel(c, 'Log Power (dB)');

end

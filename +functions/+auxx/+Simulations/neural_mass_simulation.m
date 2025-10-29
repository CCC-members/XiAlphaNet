function [cross_spectrum_avg,scalp_cross_spectrum_avg] = neural_mass_simulation(plot_index,json_path)



import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*
import functions.auxx.OptimizedOperations.*
Cortex = load("templates/Cortex.mat");

tic
% Path to the JSON file with model result metadata
[dataset_dir, ~, ~] = fileparts(json_path);
dataset = jsondecode(fileread(json_path));
dataset.Location = dataset_dir;
parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));

% Simulation parameters
N = 360;                      % Number of ROIs
tau = 2;
omega1 = 4 * pi * 0;          % 0 Hz oscillator
rng('shuffle'); % Seed with current time
freq0 = 8+5*rand(1)
omega2 = 4 * pi * freq0;         % 10 Hz oscillator
zeta1 = 0.00001;  % 0.000001 % 0.3
zeta2 = 0.00000001;
T = 1;                        % 1000 miliseconds
dt = 0.001;                   % 1 ms
n_steps = round(T / dt);
noise_std = 0.1;            % Noise standard deviation (cortical)

% Connectivity and delays
rng(0);
C = parameters.Compact_Model.C;
D = parameters.Compact_Model.D;
D_steps = round(D / dt);
R = parameters.Model.R;
K = parameters.Compact_Model.K;  % Lead field matrix

% Number of ensemble simulations
M = 20;

% Region labels
region_labels = {Cortex.Atlas(14).Scouts.Region};

% Create masks
mask1 = 0.1*ones(N, 1);  % global activation

% Create mask2: gradient increasing towards posterior regions
mask2 = zeros(N, 1);
for i = 1:N
    label = region_labels{i};
    if contains(label, 'O')  % Occipital
        mask2(i) = 0.20;
    elseif contains(label, 'P') && ~contains(label, 'F')  % Parietal
        mask2(i) = 0.17;
    elseif contains(label, 'T') && ~contains(label, 'F')  % Temporal
        mask2(i) = 0.15;
    elseif contains(label, 'C')  % Central
        mask2(i) = 0.05;
    else  % Frontal and others
        mask2(i) = 0.0;
    end
end
mask2 = 0.02 * mask2;

% Initialize ensemble storage
u_total_ensemble = zeros(M, N, n_steps);
disp('-->> Neural Mass Simulatation Initiated')

% Simulation loop
parfor m = 1:M
    u1 = zeros(N, n_steps);
    v1 = zeros(N, n_steps);
    u2 = zeros(N, n_steps);
    v2 = zeros(N, n_steps);
    
    % Seed initial condition with masks
    u1(:,1) = mask1;
    u2(:,1) = mask2;

    % Noise input per oscillator
    I1 = noise_std * randn(N, n_steps);
    I2 = noise_std * randn(N, n_steps);

    for t = 2:n_steps
        for i = 1:N
            input_sum = 0;
            for j = 1:N
                delay_idx = t - D_steps(i, j);
                if delay_idx < 1
                    delayed_activity = 0;
                else
                    delayed_activity = u1(j, delay_idx) + u2(j, delay_idx); 
                end
                input_sum = input_sum + C(i, j) * (1 / (1 + exp(-delayed_activity)));
            end

            % Oscillator 1 (0 Hz, global)
            du1 = v1(i, t-1);
            dv1 = -2 * zeta1 * omega1 * v1(i, t-1) - omega1^2 * u1(i, t-1) + input_sum + I1(i, t);
            u1(i, t) = u1(i, t-1) + dt * du1 / tau;
            v1(i, t) = v1(i, t-1) + dt * dv1 / tau;

            % Oscillator 2 (10 Hz, posterior gradient)
            du2 = v2(i, t-1);
            dv2 = -2 * zeta2 * omega2 * v2(i, t-1) - omega2^2 * u2(i, t-1) + input_sum + I2(i, t);
            u2(i, t) = u2(i, t-1) + dt * du2 / tau;
            v2(i, t) = v2(i, t-1) + dt * dv2 / tau;
        end
    end

    % Sum outputs
    u_total_ensemble(m, :, :) = u1 + u2;
end

disp('-->> Neural Mass Simulation Completed');

% FFT frequencies
freqs_full = (0:(n_steps/2)) / (n_steps*dt);

% Choose 47 evenly spaced frequency indices between 0 and 20 Hz
Nfreq = 47;
target_freqs = linspace(0, 20, Nfreq);
fft_freq_indices = zeros(1, Nfreq);
for k = 1:Nfreq
    [~, idx] = min(abs(freqs_full - target_freqs(k)));
    fft_freq_indices(k) = idx;
end
bin_centers = freqs_full(fft_freq_indices);

% Initialize averages
power_spectrum_avg = zeros(N, Nfreq);
cross_spectrum_avg = zeros(N, N, Nfreq);

% Cortex: compute power and cross-spectra
for m = 1:M
    u_total = squeeze(u_total_ensemble(m, :, :));    % [N x n_steps]
    U_fft = fft(u_total, [], 2);
    U_fft = U_fft(:, 1:(n_steps/2+1));  % positive freqs

    % Select desired frequency bins
    U_fft_selected = U_fft(:, fft_freq_indices);

    % Power spectrum
    P_selected = abs(U_fft_selected).^2 / n_steps;
    power_spectrum_avg = power_spectrum_avg + P_selected;

    % Cross-spectrum (vectorized)
    for k = 1:Nfreq
        S_k = (U_fft_selected(:, k) * U_fft_selected(:, k)') / n_steps;
        cross_spectrum_avg(:, :, k) = cross_spectrum_avg(:, :, k) + S_k;
    end
end
power_spectrum_avg = power_spectrum_avg / M;
cross_spectrum_avg = cross_spectrum_avg / M;

% Scalp: compute power and cross-spectra
n_sensors = size(K, 1);
measurement_noise_std = 0.001;

scalp_power_spectrum_avg = zeros(n_sensors, Nfreq);
scalp_cross_spectrum_avg = zeros(n_sensors, n_sensors, Nfreq);

for m = 1:M
    u_total = squeeze(u_total_ensemble(m, :, :));     % [N x n_steps]
    scalp_u_total = K * u_total;                      % [sensors x n_steps]
    scalp_u_total = scalp_u_total + measurement_noise_std * randn(n_sensors, n_steps);

    U_fft_scalp = fft(scalp_u_total, [], 2);
    U_fft_scalp = U_fft_scalp(:, 1:(n_steps/2+1));    % positive freqs

    % Select desired frequency bins
    U_fft_selected_scalp = U_fft_scalp(:, fft_freq_indices);

    % Power spectrum
    P_selected_scalp = abs(U_fft_selected_scalp).^2 / n_steps;
    scalp_power_spectrum_avg = scalp_power_spectrum_avg + P_selected_scalp;

    % Cross-spectrum (vectorized)
    for k = 1:Nfreq
        S_k_scalp = (U_fft_selected_scalp(:, k) * U_fft_selected_scalp(:, k)') / n_steps;
        scalp_cross_spectrum_avg(:, :, k) = scalp_cross_spectrum_avg(:, :, k) + S_k_scalp;
    end
end
scalp_power_spectrum_avg = scalp_power_spectrum_avg / M;
scalp_cross_spectrum_avg = scalp_cross_spectrum_avg / M;

disp('-->> Cross-spectrum computation completed.');



% === PLOTS ===
if plot_index == 1
    % Plot log10 power spectrum (cortex)
    figure('Color', 'w');
    hold on;
    colors = lines(N);
    for i = 1:N
        plot(bin_centers, 10*log10(power_spectrum_avg(i, :)), 'Color', colors(i, :), 'LineWidth', 1.5);
    end
    xlim([0 20]);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('Log10 Power Spectrum (Cortex, Ensemble Averaged)');
    grid on;
    
    % Plot cross-spectrum magnitude at 10 Hz (cortex)
    [~, idx10Hz] = min(abs(bin_centers - 10));
    cross_spectrum_mag = abs(cross_spectrum_avg(:, :, idx10Hz));
    
    figure('Color', 'w');
    imagesc(cross_spectrum_mag);
    colorbar;
    xlabel('ROI j');
    ylabel('ROI i');
    title(sprintf('Cortex Cross-Spectrum Magnitude at %.2f Hz', bin_centers(idx10Hz)));
    axis square;
    colormap(hot);
    
    % Plot log10 power spectrum (scalp)
    figure('Color', 'w');
    hold on;
    colors = lines(n_sensors);
    for i = 1:n_sensors
        plot(bin_centers, 10*log10(scalp_power_spectrum_avg(i, :)), 'Color', colors(i, :), 'LineWidth', 1.5);
    end
    xlim([0 20]);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('Log10 Power Spectrum (Scalp, Ensemble Averaged)');
    grid on;
    
    % Plot cross-spectrum magnitude at 10 Hz (scalp)
    scalp_cross_spectrum_mag = abs(scalp_cross_spectrum_avg(:, :, idx10Hz));
    
    figure('Color', 'w');
    imagesc(scalp_cross_spectrum_mag);
    colorbar;
    xlabel('Sensor j');
    ylabel('Sensor i');
    title(sprintf('Scalp Cross-Spectrum Magnitude at %.2f Hz', bin_centers(idx10Hz)));
    axis square;
    colormap(hot);
    
    toc;
    
    J = R'*power_spectrum_avg(: ,idx10Hz);
    guide.Visualization.esi_plot_single
    title('Alpha-Process')
    
    J = R'*power_spectrum_avg(:, 1);
    guide.Visualization.esi_plot_single
    title('Xi-Process')
end
end
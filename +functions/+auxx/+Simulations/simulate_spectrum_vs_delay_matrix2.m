% Hybrid Cross-Spectrum with 3x3 Visualization

clear; clc;

% --- Parameters ---
Nv = 300;
positions = linspace(0, 100, Nv);
ref_voxel = 150;
ref_pos = positions(ref_voxel);
distances = positions - ref_pos;

freqs = linspace(1, 40, 200);
Nf = length(freqs);
epsilon = 0.01;
I = eye(Nv);
distance_thresh = 2;  % mm: use diagonal approximation for closer than this

% --- Connectivity ---
C_real = exp(-abs(positions' - positions) / 15);
C_real = C_real + 0.05 * rand(Nv);
C_real(1:Nv+1:end) = 0;

% --- Delays ---
D_real = abs(positions' - positions) / 10 + 0.1 * rand(Nv);
D_real(1:Nv+1:end) = 0;

% --- Spectral Components ---
xi_power = 10 ./ (1 + (freqs(:) / 7).^2).^3.5;
xi_power = repmat(xi_power, 1, Nv);

alpha_power = 10 ./ (1 + ((freqs(:) - 10)/5).^2).^2;
alpha_power = repmat(alpha_power, 1, Nv);

total_power = xi_power + alpha_power;

% --- Hybrid Cross-Spectrum Computation ---
cross_total = zeros(Nf, Nv);
cross_xi    = zeros(Nf, Nv);
cross_alpha = zeros(Nf, Nv);

for f_idx = 1:Nf
    omega = 2 * pi * freqs(f_idx);
    A_omega = zeros(Nv);
    for k = 1:Nv
        if k ~= ref_voxel
            A_omega(ref_voxel, k) = C_real(ref_voxel, k) * exp(-1i * omega * D_real(ref_voxel, k));
        end
    end
    G_omega = inv((1 + epsilon) * I - A_omega);
    g_i = G_omega(ref_voxel, :);

    for j = 1:Nv
        d = abs(positions(j) - ref_pos);
        if d < distance_thresh
            % Diagonal approximation
            Gval = abs(G_omega(ref_voxel, j))^2;
            cross_total(f_idx, j) = Gval * total_power(f_idx, j);
            cross_xi(f_idx, j)    = Gval * xi_power(f_idx, j);
            cross_alpha(f_idx, j) = Gval * alpha_power(f_idx, j);
        else
            % Full cross-spectrum
            g_j = G_omega(j, :);
            cross_total(f_idx, j) = abs(sum(g_i .* conj(g_j) .* total_power(f_idx, :)));
            cross_xi(f_idx, j)    = abs(sum(g_i .* conj(g_j) .* xi_power(f_idx, :)));
            cross_alpha(f_idx, j) = abs(sum(g_i .* conj(g_j) .* alpha_power(f_idx, :)));
        end
    end
end

% Convert to log scale
log_cross_total = 10 * log10(cross_total + 1e-10);
log_cross_xi    = 10 * log10(cross_xi + 1e-10);
log_cross_alpha = 10 * log10(cross_alpha + 1e-10);

% Grid
[DIST, FREQ] = meshgrid(distances, freqs);

% --- 3x3 Plot Layout ---
figure('Position', [100, 100, 1400, 1000]);
sgtitle('Hybrid Cross-Spectrum Components (Diagonal for Short Range)', 'FontSize', 16);

% Row 1: Total
subplot(3,3,1); surf(DIST, FREQ, log_cross_total, 'EdgeColor', 'none');
view([45 30]); title('Total - Side View');
xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power'); colormap parula;

subplot(3,3,2); surf(DIST, FREQ, log_cross_total, 'EdgeColor', 'none');
view([120 45]); title('Total - Oblique View');
xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');

subplot(3,3,3); imagesc(freqs, distances, log_cross_total'); axis xy;
xlabel('Frequency'); ylabel('Distance'); title('Total - Heatmap');
colorbar;

% Row 2: Xi
subplot(3,3,4); surf(DIST, FREQ, log_cross_xi, 'EdgeColor', 'none');
view([45 30]); title('\xi - Side View');
xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');

subplot(3,3,5); surf(DIST, FREQ, log_cross_xi, 'EdgeColor', 'none');
view([120 45]); title('\xi - Oblique View');
xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');

subplot(3,3,6); imagesc(freqs, distances, log_cross_xi'); axis xy;
xlabel('Frequency'); ylabel('Distance'); title('\xi - Heatmap');
colorbar;

% Row 3: Alpha
subplot(3,3,7); surf(DIST, FREQ, log_cross_alpha, 'EdgeColor', 'none');
view([45 30]); title('\alpha - Side View');
xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');

subplot(3,3,8); surf(DIST, FREQ, log_cross_alpha, 'EdgeColor', 'none');
view([120 45]); title('\alpha - Oblique View');
xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');

subplot(3,3,9); imagesc(freqs, distances, log_cross_alpha'); axis xy;
xlabel('Frequency'); ylabel('Distance'); title('\alpha - Heatmap');
colorbar;

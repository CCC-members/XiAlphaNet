function hybrid_cross_spectrum_kernel(C, D)
% HYBRID_CROSS_SPECTRUM_KERNEL
%   Compute and visualize the hybrid cross-spectrum using a precomputed kernel
%   from connectivity (C) and delay (D) matrices.

    Nv = size(C, 1);
    positions = linspace(0, 100, Nv);
    ref_voxel = round(Nv/2); 
    ref_pos = positions(ref_voxel);
    distances = positions - ref_pos;

    freqs = linspace(1, 40, 200);
    Nf = length(freqs);
    epsilon = 0.01;
    I = eye(Nv);
    distance_thresh = 2;  % mm

    % --- Spectral Components ---
    xi_power = 10 ./ (1 + (freqs(:) / 7).^2).^3.5;
    xi_power = repmat(xi_power, 1, Nv);

    alpha_power = 10 ./ (1 + ((freqs(:) - 10)/5).^2).^2;
    alpha_power = repmat(alpha_power, 1, Nv);

    total_power = xi_power + alpha_power;

    % --- Initialize Results ---
    cross_total = zeros(Nf, Nv);
    cross_xi    = zeros(Nf, Nv);
    cross_alpha = zeros(Nf, Nv);

    % --- Main Loop over Frequency ---
    for f_idx = 1:Nf
        omega = 2 * pi * freqs(f_idx);

        % === Build complex kernel ===
        C_kernel = C .* exp(-1i * omega * D);
        C_kernel(1:Nv+1:end) = 0;  % force zero diagonal

        % === Compute Green's function ===
        G_omega = inv((1 + epsilon) * I - C_kernel);
        g_i = G_omega(ref_voxel, :);  % Row for source

        for j = 1:Nv
            d = abs(positions(j) - ref_pos);
            if d < distance_thresh
                Gval = abs(G_omega(ref_voxel, j))^2;
                cross_total(f_idx, j) = Gval * total_power(f_idx, j);
                cross_xi(f_idx, j)    = Gval * xi_power(f_idx, j);
                cross_alpha(f_idx, j) = Gval * alpha_power(f_idx, j);
            else
                g_j = G_omega(j, :);
                cross_total(f_idx, j) = abs(sum(g_i .* conj(g_j) .* total_power(f_idx, :)));
                cross_xi(f_idx, j)    = abs(sum(g_i .* conj(g_j) .* xi_power(f_idx, :)));
                cross_alpha(f_idx, j) = abs(sum(g_i .* conj(g_j) .* alpha_power(f_idx, :)));
            end
        end
    end

    % --- Log Power ---
    log_cross_total = 10 * log10(cross_total + 1e-10);
    log_cross_xi    = 10 * log10(cross_xi + 1e-10);
    log_cross_alpha = 10 * log10(cross_alpha + 1e-10);

    % --- Plot ---
    [DIST, FREQ] = meshgrid(distances, freqs);
    figure('Position', [100, 100, 1400, 1000]);
    sgtitle('Hybrid Cross-Spectrum using Complex Kernel', 'FontSize', 16);

    % Total
    subplot(3,3,1); surf(DIST, FREQ, log_cross_total, 'EdgeColor', 'none');
    view([45 30]); title('Total - Side View');
    xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power'); colormap parula;
    subplot(3,3,2); surf(DIST, FREQ, log_cross_total, 'EdgeColor', 'none');
    view([120 45]); title('Total - Oblique View');
    xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');
    subplot(3,3,3); imagesc(freqs, distances, log_cross_total'); axis xy;
    xlabel('Frequency'); ylabel('Distance'); title('Total - Heatmap'); colorbar;

    % Xi
    subplot(3,3,4); surf(DIST, FREQ, log_cross_xi, 'EdgeColor', 'none');
    view([45 30]); title('\xi - Side View');
    xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');
    subplot(3,3,5); surf(DIST, FREQ, log_cross_xi, 'EdgeColor', 'none');
    view([120 45]); title('\xi - Oblique View');
    xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');
    subplot(3,3,6); imagesc(freqs, distances, log_cross_xi'); axis xy;
    xlabel('Frequency'); ylabel('Distance'); title('\xi - Heatmap'); colorbar;

    % Alpha
    subplot(3,3,7); surf(DIST, FREQ, log_cross_alpha, 'EdgeColor', 'none');
    view([45 30]); title('\alpha - Side View');
    xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');
    subplot(3,3,8); surf(DIST, FREQ, log_cross_alpha, 'EdgeColor', 'none');
    view([120 45]); title('\alpha - Oblique View');
    xlabel('Distance'); ylabel('Frequency'); zlabel('Log Power');
    subplot(3,3,9); imagesc(freqs, distances, log_cross_alpha'); axis xy;
    xlabel('Frequency'); ylabel('Distance'); title('\alpha - Heatmap'); colorbar;
end

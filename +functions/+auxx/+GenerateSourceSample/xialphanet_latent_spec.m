function [Sjj] = xialphanet_latent_spec(Svv, T, useParallel)
% XIALPHANET_LATENT_SPEC_4BANDS
% --------------------------------------------------------------
% Compute latent source spectral distribution (diagonal only)
% using Xi-AlphaNET inverse operator with eLORETA-like
% iterative weighting + GCV scheme.
% Here: 4 frequency bands, each with its own W and gamma.
%
% Inputs:
%   Svv         : (Ne × Ne × Nw) sensor cross-spectra
%   T           : (Ne × Nv × Nw) forward operator per frequency
%   useParallel : logical, true=parfor
%
% Output:
%   Sjj : (Nv × Nw) source power spectra
% --------------------------------------------------------------

    [Ne, ~, Nw] = size(Svv);
    [~, Nv, Nk] = size(T);
    if Nk ~= Nw
        error('T must be (Ne × Nv × Nw)');
    end

    Sjj = zeros(Nv, Nw);

    % === Step 1. Define 4 bins ===
    edges = round(linspace(1, Nw+1, 5));  % 4 bins
    bins = cell(1,4);
    for b = 1:4
        bins{b} = edges(b):(edges(b+1)-1);
    end

    % === Step 2. For each bin compute W_b and gamma_b ===
    for b = 1:4
        idxs = bins{b};
        mid_idx = round(mean(idxs));
        Fop_mid = T(:,:,mid_idx);

        % --- Compute eLORETA weights for this band ---
        regu0 = 0.05;
        [~, Wb] = mkiter_xialphanet(Fop_mid, regu0);

        % --- Average cross-spectra in the band ---
        Svv_avg = mean(Svv(:,:,idxs), 3);

        % --- Compute gamma_b via GCV ---
        gamma_b = gcv_find_gamma(Svv_avg, Fop_mid, Wb);

        % --- Fixed operator for this band ---
        M = (Fop_mid*Wb*Fop_mid' + gamma_b*eye(Ne)) \ eye(Ne);
        UF_mid = M * Fop_mid * Wb;   % (Ne × Nv)

        % === Step 3. Apply to all frequencies in bin ===
        if useParallel
            Sj_local = zeros(Nv, numel(idxs));
            parfor ii = 1:numel(idxs)
                iw = idxs(ii);
                UF = M * T(:,:,iw) * Wb;
                Sj_local(:,ii) = abs(sum((UF' * Svv(:,:,iw)) .* UF', 2));
            end
            Sjj(:,idxs) = Sj_local;
        else
            for ii = 1:numel(idxs)
                iw = idxs(ii);
                UF = M * T(:,:,iw) * Wb;
                Sjj(:,iw) = abs(sum((UF' * Svv(:,:,iw)) .* UF', 2));
            end
        end
    end
end


%% === Iterative eLORETA weight update ===
function [M, Wout] = mkiter_xialphanet(Fop, regu)
    [Ne, Nv] = size(Fop);
    Wout = eye(Nv);
    maxIter = 20; tol = 1e-6;

    for kk = 1:maxIter
        Winv = Wout \ eye(Nv);
        kwinvkt = Fop * Winv * Fop';
        alpha   = regu * trace(kwinvkt) / Ne;

        M = (kwinvkt + alpha*eye(Ne)) \ eye(Ne);

        % Vectorized voxel weight update (correct eLORETA)
        d = sum((Fop' * M) .* Fop', 2);   % diag(F' M F)
        Wnew = 1 ./ sqrt(real(d));        % <-- FIXED

        % Update diagonal matrix W
        Wold = diag(Wout);
        Wout = diag(Wnew);

        % Convergence check
        if norm(Wnew - Wold)/norm(Wold) < tol
            break;
        end
    end
end


%% === GCV for one band (fixed W) ===
function gamma = gcv_find_gamma(Svv, Fop, W)
    Ne = size(Fop,1);
    gamma_grid = 0.01:0.01:0.5;
    gcv = zeros(numel(gamma_grid),1);

    for k = 1:numel(gamma_grid)
        gamma = gamma_grid(k);
        M = (Fop*W*Fop' + gamma*eye(Ne)) \ eye(Ne);   % <-- FIXED: use W not W^-1
        R = eye(Ne) - Fop * W * (Fop' * M);
        num = trace(R * Svv * R');
        den = (trace(R))^2;
        gcv(k) = num / den;
    end

    [~, idx] = min(gcv);
    gamma = gamma_grid(idx);
end

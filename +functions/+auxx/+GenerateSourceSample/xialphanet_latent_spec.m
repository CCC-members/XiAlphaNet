function [Sjj] = xialphanet_latent_spec(Svv, T, useParallel)
% XIALPHANET_LATENT_SPEC
% --------------------------------------------------------------
% Compute latent source spectral distribution (diagonal only)
% using the Xi-AlphaNET inverse operator in an eLORETA-like
% iterative weighting + GCV scheme (voxel-wise, no big matrices).
%
%   Inputs:
%       Svv         : (Ne × Ne × Nw) sensor cross-spectra
%       T           : (Ne × Nv × Nw) forward operator (freq-dependent)
%       useParallel : logical flag, true = parfor
%
%   Output:
%       Sjj : (Nv × Nw) latent source spectra (diagonal only)
% --------------------------------------------------------------

    [Ne, ~, Nw] = size(Svv);
    [~, Nv, Nk] = size(T);
    if Nk ~= Nw
        error('T must be (Ne × Nv × Nw)');
    end

    Sjj = zeros(Nv, Nw);

    if useParallel
        parfor iw = 1:Nw
            iw
            Sjj(:,iw) = gcv_latent_diag(Svv(:,:,iw), T(:,:,iw));
        end
    else
        for iw = 1:Nw
            Sjj(:,iw) = gcv_latent_diag(Svv(:,:,iw), T(:,:,iw));
        end
    end
end


%% ===== GCV with voxel-wise diagonal computation =====
function Sdiag = gcv_latent_diag(Svv, Fop)
    Ne = size(Fop,1);
    Nv = size(Fop,2);

    gamma_grid = 0.01:0.01:0.5;
    gcv = zeros(numel(gamma_grid),1);

    % --- Step 1. Grid search ---
    for k = 1:numel(gamma_grid)
        gamma = gamma_grid(k);
        [M, W] = mkiter_xialphanet(Fop, gamma);

        % Residual operator
        R = eye(Ne) - Fop * (W \ (Fop' * M));
        num = trace(R * Svv * R');
        den = (trace(R))^2;
        gcv(k) = num / den;
    end

    % --- Step 2. Choose optimal gamma ---
    [~, idx] = min(gcv);
    gamma = gamma_grid(idx);

    % --- Step 3. Final inverse and voxel spectra ---
    [M, W] = mkiter_xialphanet(Fop, gamma);

    Sdiag = zeros(Nv,1);
    for i = 1:Nv
        ui = M * Fop(:,i);                         % (Ne × 1)
        Sdiag(i) = (1/W(i,i)^2) * real(ui' * Svv * ui);
    end
end


%% ===== Iterative eLORETA weight update =====
function [M, Wout] = mkiter_xialphanet(Fop, regu)
    [Ne, Nv] = size(Fop);
    Wout = eye(Nv);
    maxIter = 20; tol = 1e-6;

    for kk = 1:maxIter
        Winv = Wout \ eye(Nv);
        kwinvkt = Fop * Winv * Fop';
        alpha   = regu * trace(kwinvkt) / Ne;

        M = (kwinvkt + alpha*eye(Ne)) \ eye(Ne);

        Wold = diag(Wout);
        for i = 1:Nv
            li = Fop(:,i);
            Wout(i,i) = sqrt(real(li' * M * li));
        end

        if norm(diag(Wout)-Wold)/norm(Wold) < tol
            break;
        end
    end
end

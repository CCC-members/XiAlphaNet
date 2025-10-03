function [Sjj] = eloreta_cross(Svv, K, useParallel)
% XIALPHANET_SPECTRA Compute Xi-AlphaNET source spectra (diagonal only)
% using frequency-dependent forward operators and eLORETA-like inversion
%
% Inputs:
%   Svv        : (Ne × Ne × Nw) sensor cross-spectra
%   K          : (Ne × Nv × Nw) frequency-dependent forward operator
%   useParallel: logical flag, true = use parfor, false = normal for
%
% Output:
%   Sjj : (Nv × Nw) source power spectra

    [Ne, ~, Nw] = size(Svv);
    [~, Nv, Nk] = size(K);

    if Nk ~= Nw
        error('K must be frequency-dependent: size(K,3) = Nw.');
    end

    % --- Allocate directly ---
    Sjj = zeros(Nv, Nw);

    % --- Choose parallel or serial loop ---
    if useParallel
        parfor iw = 1:Nw
            iw
            Sjj(:,iw) = gcv_xialphanet(Svv(:,:,iw), K(:,:,iw));
        end
    else
        for iw = 1:Nw
            Sjj(:,iw) = gcv_xialphanet(Svv(:,:,iw), K(:,:,iw));
        end
    end
end


%% ===== Helper: Xi-AlphaNET GCV =====
function Sdiag = gcv_xialphanet(Svv, Fop)
    p   = size(Fop,1);
    Ip  = eye(p);

    % Parameters for grid search
    gamma_grid = 0.01:0.01:0.5;
    gcv = zeros(numel(gamma_grid),1);

    % Grid search over gamma
    for k = 1:numel(gamma_grid)
        gamma = gamma_grid(k);
        Tjv   = mkfilt_xialphanet(Fop, gamma);   % (Nv × Ne)
        Txiv  = Ip - Fop*Tjv;                    % (Ne × Ne)
        num   = sum(abs(diag(Txiv*Svv*Txiv')));
        den   = sum(abs(diag(Txiv)))^2;
        gcv(k) = num / den;
    end

    % Select optimal gamma
    [~, idx] = min(gcv);
    gamma    = gamma_grid(idx);

    % Final filter and spectra
    Tjv   = mkfilt_xialphanet(Fop, gamma);
    Sdiag = real(diag(Tjv * Svv * Tjv'));
end


%% ===== Xi-AlphaNET Filter (adapted eLORETA) =====
function A = mkfilt_xialphanet(Fop, regu, W)
    if nargin < 2, regu = 0.05; end
    [Ne, Nv] = size(Fop);

    if nargin < 3
        W = eye(Nv);
    end

    u0 = eye(Ne);
    Wout = W;
    kont = false; kk = 0;

    while ~kont
        kk = kk + 1;

        % === Mixed precision inversion in voxel space ===
        Wtmp = single(Wout + trace(Wout)/(Nv*1e6));
        Winv = double(Wtmp \ eye(Nv,'single'));

        % --- Build system matrix (sensor space, smaller) ---
        kwinvkt = Fop * Winv * Fop';
        alpha   = regu * trace(kwinvkt) / Ne;

        % === Mixed precision inversion in sensor space ===
        Mtmp = single(kwinvkt + alpha*u0);
        M    = double(Mtmp \ eye(Ne,'single'));

        % --- Update voxel weights ---
        Wold = Wout;
        for i = 1:Nv
            li = Fop(:,i);
            Mb = li' * M * li;
            Wout(i,i) = sqrt(real(Mb) + trace(Mb)/1e10);
        end

        % --- Convergence check ---
        reldef = norm(Wout(:)-Wold(:)) / norm(Wold(:));
        if kk > 20 || reldef < 1e-6
            kont = true;
        end
    end

    % --- Final spatial filter ---
    A = Wout \ (Fop' * M);   % (Nv × Ne)
end

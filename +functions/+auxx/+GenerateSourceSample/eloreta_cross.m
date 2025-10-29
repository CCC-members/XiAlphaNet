function [Sjj] = eloreta_cross(Svv, K, useParallel, R)
% ELORETA_CROSS  Compute Xi-AlphaNET source spectra (diagonal only)
% using frequency-dependent forward operators and eLORETA-like inversion.
%
% Inputs:
%   Svv        : (Ne × Ne × Nw) sensor cross-spectra
%   K          : (Ne × Nv × Nw) frequency-dependent forward operator
%   useParallel: logical flag, true = use parfor, false = normal for
%   R          : (Nr × Nv) ROI?voxel projection (optional, accelerates inversion)
%
% Output:
%   Sjj : (Nv × Nw) voxel-wise source power spectra

[Ne, ~, Nw] = size(Svv);
[~, Nv, Nk] = size(K);
if Nk ~= Nw
    error('K must be frequency-dependent: size(K,3) = Nw.');
end

Sjj = zeros(Nv, Nw);

% -------------------------------------------------------------------------
%                     PARALLEL EXECUTION WITH TWO BATCHES
% -------------------------------------------------------------------------
if useParallel
    midN    = ceil(Nw/2);
    batches = {1:midN, midN+1:Nw};

    for b = 1:2
        idx = batches{b};
        nb  = numel(idx);
        fprintf('\n--- Processing batch %d/%d (%d frequencies) ---\n', b, numel(batches), nb);

        % Temporary container (parfor-safe)
        Sbatch = cell(1, nb);

        parfor iw = 1:nb
            widx = idx(iw);
            Sbatch{iw} = gcv_xialphanet(Svv(:,:,widx), K(:,:,widx), R);
        end

        % Gather back into main matrix
        for iw = 1:nb
            widx = idx(iw);
            Sjj(:, widx) = Sbatch{iw};
        end
    end

else
    for iw = 1:Nw
        Sjj(:, iw) = gcv_xialphanet(Svv(:,:,iw), K(:,:,iw), R);
    end
end

end % ======================= end of main function =========================



%% ========================================================================
%                     Helper: Xi-AlphaNET GCV  (fixed ?)
% ========================================================================
function Sdiag = gcv_xialphanet(Svv, Fop, R)
gamma = 0.05;  % theoretical optimal regularization
Tjv   = mkfilt_xialphanet(Fop, gamma, [], R);
Sdiag = real(diag(Tjv * Svv * Tjv'));
end



%% ========================================================================
%                    Helper: Xi-AlphaNET Filter (adapted eLORETA)
% ========================================================================
function A = mkfilt_xialphanet(Fop, regu, W, R)
if nargin < 2 || isempty(regu), regu = 0.05; end
[Ne, Nv] = size(Fop);

if nargin < 3 || isempty(W)
    W = eye(Nv);
end

u0   = eye(Ne);
Wout = eye(Nv);
kont = false; kk = 0;

while ~kont
    kk = kk + 1;

    % === Accelerated inversion in ROI space using R (ROI × voxel) ===
    if exist('R','var') && ~isempty(R)
        % Project W to ROI space (Nr × Nr)
        W_R = R * Wout * R';
        % Regularize for stability
        W_R = W_R + trace(W_R)/(size(W_R,1)*1e6);
        % Invert in ROI space
        W_Rinv = W_R \ eye(size(W_R));
        % Project back to voxel space
        Winv = R' * W_Rinv * R;
    else
        % Standard voxel-space inversion
        Wtmp = (Wout + trace(Wout)/(Nv*1e6));
        Winv = (Wtmp \ eye(Nv));
    end

    % --- Build system matrix in sensor space ---
    kwinvkt = Fop * Winv * Fop';
    alpha   = real(regu * trace(kwinvkt) / Ne);

    % --- Invert in sensor space (Ne × Ne) ---
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
    if kk > 20 || reldef < 1e-2
        kont = true;
    end
end

% --- Final spatial filter ---
A = Wout \ (Fop' * M);   % (Nv × Ne)
end

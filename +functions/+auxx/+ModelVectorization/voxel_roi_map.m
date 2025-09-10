function R = voxel_roi_map(Cortex)
% VOXEL_ROI_MAP  Build ROI?vertex operator with orthonormal rows (R*R' = I)
%   R = voxel_roi_map(Cortex)
% Inputs:
%   Cortex.Vertices [N x 3], Cortex.Faces [F x 3]
%   Cortex.Atlas(Cortex.iAtlas).Scouts(i).Vertices : indices of vertices in ROI i
% Output:
%   R [m x N] sparse, with R*R'  I_m (area-weighted, orthonormal rows)
%
% Method:
%   R = D_A^{-1/2} * B * W^{1/2}
%   B(i,j)=1 if vertex j belongs to ROI i; W=diag(vertex areas); a_i=(B*w)_i

    Atlas = Cortex.Atlas(Cortex.iAtlas);
    N = size(Cortex.Vertices, 1);
    m = numel(Atlas.Scouts);

    % --- Build sparse incidence B (m x N): row i marks vertices in ROI i
    % Preallocate index buffers for speed
    counts = arrayfun(@(s) numel(s.Vertices), Atlas.Scouts(:));
    nnzB = sum(counts);
    I = zeros(1, nnzB);
    J = zeros(1, nnzB);
    k = 1;
    for i = 1:m
        idx = Atlas.Scouts(i).Vertices(:)';   % row vector of vertex indices
        len = numel(idx);
        I(k:k+len-1) = i;
        J(k:k+len-1) = idx;
        k = k + len;
    end
    B = sparse(I, J, 1, m, N);

    % --- Per-vertex areas (Voronoi by 1/3 triangle split)
    V = Cortex.Vertices;
    if isfield(Cortex, 'Faces') && ~isempty(Cortex.Faces)
        F = Cortex.Faces;
        v1 = V(F(:,1),:); v2 = V(F(:,2),:); v3 = V(F(:,3),:);
        Atri = 0.5 * sqrt(sum(cross(v2 - v1, v3 - v1, 2).^2, 2));   % triangle areas
        w = accumarray(F(:), repmat(Atri/3, 3, 1), [N,1], @sum, 0); % vertex areas
        w(w <= 0) = eps;
    else
        % Fallback: uniform areas if Faces are unavailable
        w = ones(N,1);
    end

    % --- ROI areas and row normalization
    a = B * w;                            % m x 1 total area per ROI
    a(a <= 0) = eps;                      % guard empty ROIs

    % R = D_A^{-1/2} * B * W^{1/2}
    R = spdiags(1 ./ sqrt(a), 0, m, m) * B * spdiags(sqrt(w), 0, N, N);

    % (Optional) quick sanity check  comment out in production
    % fprintf('||R*R'' - I||_F = %.2e\n', norm(full(R*R') - eye(m), 'fro'));
end

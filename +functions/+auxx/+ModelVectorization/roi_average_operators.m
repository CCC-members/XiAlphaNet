function [R_fwd, R_inv] = roi_average_operators(Cortex, CortexiAtlas)
% ROI_AVERAGE_OPERATORS  Build forward/inverse ROI operators without smoothing.
% 
% Inputs:
%   Cortex.Atlas(CortexiAtlas).Scouts(i).Vertices : indices of vertices in ROI i
%   N = #vertices, m = #ROIs
%
% Outputs:
%   R_fwd [m x N]  : ROI averaging operator (row i averages vertices in ROI i)
%   R_inv [N x m]  : ROI expansion operator (column i assigns ROI mean to vertices)

    Atlas = Cortex.Atlas(CortexiAtlas);
    N = size(Cortex.Vertices,1);
    m = numel(Atlas.Scouts);

    % --- Forward operator: average within ROI ---
    I = []; J = []; S = [];
    for i = 1:m
        verts = Atlas.Scouts(i).Vertices(:);
        nverts = numel(verts);
        I = [I; repmat(i,nverts,1)];
        J = [J; verts];
        S = [S; ones(nverts,1)/nverts];   % equal weights 1/|V_i|
    end
    R_fwd = sparse(I,J,S,m,N);

    % --- Inverse operator: broadcast ROI values to vertices ---
    I = []; J = []; S = [];
    for i = 1:m
        verts = Atlas.Scouts(i).Vertices(:);
        nverts = numel(verts);
        I = [I; verts];
        J = [J; repmat(i,nverts,1)];
        S = [S; ones(nverts,1)];
    end
    R_inv = sparse(I,J,S,N,m);
end

function [E] = roi_voxel_map()
 % Load the anatomical atlas data (adjust the path if needed)
    load('Data/Atlas_Anatomical/tess_cortex_mid_high_8000V_fix.mat');
    HCP_MMP1 = Atlas(iAtlas);  % Extract the atlas (assuming iAtlas is defined)
    Nv = length(Curvature);  % Number of voxels in the brain
    Nr = 360;  % Assuming 360 ROIs based on your anatomical atlas

    % Initialize the voxel-to-ROI projection matrix E2
    E = zeros(Nv, Nr);
    
    % For each ROI, map its vertices (voxels) and assign normalized values in E2
    for j = 1:Nr
        vertices_j = HCP_MMP1.Scouts(j).Vertices;  % Vertices (voxels) for ROI j
        n_j = length(vertices_j);  % Number of voxels in this ROI
        
        % Assign values to each voxel in the ROI
        % Normalize by sqrt(n_j) to ensure that E2' * E2 = I
        E(vertices_j, j) = 1 / sqrt(n_j);
    end
end

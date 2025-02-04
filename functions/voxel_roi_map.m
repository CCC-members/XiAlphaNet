function [R] = voxel_roi_map
    load('Data/Atlas_Anatomical/tess_cortex_mid_high_8000V_fix.mat');
    HCP_MMP1 = Atlas(iAtlas);
    N_voxel = length(Curvature);
    
    % Initialize R with zeros
    R = zeros(360, N_voxel);
    
    for i = 1:360
        vertices_i = HCP_MMP1.Scouts(i).Vertices;
        
        % Calculate the mapping values for the current ROI
        for j = 1:N_voxel
            R(i, j) = sum(vertices_i == j);
        end
    end
    
    % Normalize the rows of R
    for i = 1:size(R, 1)
        norm_i = norm(R(i, :));  % Calculate the norm of the row
        if norm_i > 0
            R(i, :) = R(i, :) / norm_i;  % Normalize the row
        end
    end
end




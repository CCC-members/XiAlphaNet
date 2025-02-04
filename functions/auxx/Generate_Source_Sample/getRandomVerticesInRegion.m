% Function to select random vertices from a specific region (either left or right hemisphere)
function vertices = getRandomVerticesInRegion(hemisphere)
    vertices = [];
    
    % Define regions for each lobe and hemisphere using the provided names
    if strcmp(hemisphere, 'left')
        % Left hemisphere regions
        occipital_regions = {'L_V1_ROI L', 'L_V2_ROI L', 'L_V3A_ROI L', 'L_V3B_ROI L', 'L_V3CD_ROI L', 'L_V3_ROI L', 'L_V4_ROI L', 'L_V7_ROI L'};
        cuneus_regions    = {'L_V6_ROI L', 'L_V6A_ROI L', 'L_V8_ROI L'};
        parietal_regions  = {'L_7m_ROI L', 'L_7PC_ROI L', 'L_7AL_ROI L', 'L_7PL_ROI L', 'L_7Pm_ROI L', 'L_7Am_ROI L', 'L_PGi_ROI L', 'L_PGp_ROI L', 'L_PGs_ROI L'};
    elseif strcmp(hemisphere, 'right')
        % Right hemisphere regions
        occipital_regions = {'R_V1_ROI R', 'R_V2_ROI R', 'R_V3A_ROI R', 'R_V3B_ROI R', 'R_V3CD_ROI R', 'R_V3_ROI R', 'R_V4_ROI R', 'R_V7_ROI R'};
        cuneus_regions    = {'R_V6_ROI R', 'R_V6A_ROI R', 'R_V8_ROI R'};
        parietal_regions  = {'R_7m_ROI R', 'R_7PC_ROI R', 'R_7AL_ROI R', 'R_7PL_ROI R', 'R_7Pm_ROI R', 'R_7Am_ROI R', 'R_PGi_ROI R', 'R_PGp_ROI R', 'R_PGs_ROI R'};
    end
    
    % Combine all regions into a single list
    all_regions = [occipital_regions, cuneus_regions, parietal_regions];
    
    % Select a random region from the list
    selected_region = all_regions{randi(length(all_regions))};
    
    % Get the vertices for the selected region (assuming getRegionVertices does this)
    vertices = getRegionVertices(selected_region);
end

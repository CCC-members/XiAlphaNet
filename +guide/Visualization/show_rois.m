% Define verified regions for each lobe in the HCP Atlas
frontalRegions = {'10d_ROI', '10pp_ROI', '10r_ROI', 'a10p_ROI', '9m_ROI', '9p_ROI', '11l_ROI'};
temporalRegions = {'TE1a_ROI', 'TE1p_ROI', 'TE2a_ROI', 'TE2p_ROI', 'STSvp_ROI', 'STSdp_ROI', 'PHT_ROI', 'PH_ROI'};
parietalRegions = {'7AL_ROI', '7PC_ROI', '7PL_ROI', '7Pm_ROI', 'PGp_ROI', 'PGs_ROI', 'PFm_ROI', 'PF_ROI'};
occipitalRegions = {'V1_ROI', 'V2_ROI', 'V3A_ROI', 'V3B_ROI', 'V4_ROI', 'V8_ROI', 'LO1_ROI', 'LO2_ROI'};

% Define hemispheres
hemispheres = {'L', 'R'}; % Left and Right hemispheres

% Initialize vector J
J = zeros(8003, 1); % Initialize the vertex assignment vector

% Assign numbers for each lobe
lobeNumbers = struct('Frontal', 1, 'Temporal', 1, 'Parietal', 1, 'Occipital', 1);

% Loop through each lobe and assign vertex numbers
lobes = {'Frontal', 'Temporal', 'Parietal', 'Occipital'};
regionsByLobe = {frontalRegions, temporalRegions, parietalRegions, occipitalRegions};

for l = 1:length(lobes)
    lobeName = lobes{l}; % Get the lobe name
    lobeRegions = regionsByLobe{l}; % Get the regions for this lobe
    lobeNumber = lobeNumbers.(lobeName); % Get the number assigned to this lobe
    
    for h = 1:length(hemispheres)
        hemisphere = hemispheres{h}; % 'L' or 'R'
        
        for r = 1:length(lobeRegions)
            % Construct region name (e.g., 'L_10d_ROI L' or 'R_10d_ROI R')
            regionName = sprintf('%s_%s %s', hemisphere, lobeRegions{r}, hemisphere);
            
            % Get vertices for this region
            try
                regionVertices = getRegionVertices(regionName); % Get vertices for the region
                % Assign the lobe number to these vertices
                J(regionVertices) = lobeNumber;
            catch
                % If region is not found, display a warning and skip
                warning('Region "%s" not found in the atlas. Skipping...', regionName);
            end
        end
    end
end

% Display or process `J` as needed
disp('Vertex assignments have been completed.');

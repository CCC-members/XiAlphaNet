function vertices = getRegionVertices(regionName)
    % getRegionVertices - Returns the vertices of the specified region
    %
    % Syntax: vertices = getRegionVertices(regionName)
    %
    % Inputs:
    %    regionName - The name of the region for which vertices are needed
    %
    % Outputs:
    %    vertices - Array of vertices corresponding to the regionName
    
    % Load the atlas structure
    import templates.*  
    loadedData = load('+templates/Cortex.mat');
    
    % Assuming the atlas is stored in Atlas(8).Scouts
    atlas = loadedData.Atlas(8).Scouts;
    
    % Initialize output
    vertices = [];
    
    % Loop through the atlas scouts to find the matching region name
    for i = 1:length(atlas)
        if strcmp(atlas(i).Label, regionName)
            vertices = atlas(i).Vertices;
            return;  % Exit the function once the region is found
        end
    end
    
    % If no matching region is found
    if isempty(vertices)
        warning('Region "%s" not found in the atlas.', regionName);
    end
end

function esi_plot(ax, J, colorLimits)
    % Include all your existing code here, adjusting plot commands to use 'ax'

    %% Adding necessary paths
    addpath("Visualization/ESI_plot/ESI_plot/app/");
    addpath("Visualization/ESI_plot/ESI_plot/data/");
    addpath("Visualization/ESI_plot/ESI_plot/functions/");
    addpath("Visualization/ESI_plot/ESI_plot/template/");

    %% Load data
    Cortex = load("tess_cortex_mid_high_8000V_fix");
    template = load("template/axes.mat");
    colorMap = load("template/mycolormap_brain_basic_conn.mat");
    hemis = []; % left or right or []
    [Cortex, iHideVert] = split_hemisphere(Cortex, hemis);
    Source = J; % Ensure J is defined in the calling script
    J(iHideVert) = [];
    sources_iv = J; 
    sources_iv = (sources_iv);
    smoothValue = 0.5;
    SurfSmoothIterations = 10;
    Vertices = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);

    %% Plotting section
    patch(ax, ...
    'Faces',Cortex.Faces, ...
    'Vertices',Vertices, ...
    'FaceVertexCData', log(1+sources_iv), ...
    'FaceColor','interp', ...
    'EdgeColor','none', ...
    'AlphaDataMapping', 'none', ...
    'EdgeColor', 'none', ...
    'EdgeAlpha', 1, ...
    'BackfaceLighting', 'lit', ...
    'AmbientStrength', 0.5, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.2, ...
    'SpecularExponent', 1, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud', ...
    'EdgeLighting', 'gouraud', ...
    'FaceAlpha', .99);
    view(ax, 0, 0); % Adjust view for 3D plots if necessary
    set(ax,'xcolor','w','ycolor','w','zcolor','w');
    colormap(ax, colorMap.cmap_a);
    axis(ax, 'equal');
    rotate3d(ax, 'on');
    axis off
  
    %% Set the color limits for comparison
    if exist('colorLimits', 'var') && ~isempty(colorLimits)
        caxis(ax, colorLimits);
    else
        caxis(ax, 'auto'); % Or manually set a default range if necessary
    end
end


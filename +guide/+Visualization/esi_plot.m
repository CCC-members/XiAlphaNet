function esi_plot(ax, J, colorLimits)
    % Ensure the colormap is applied to the plot by directly assigning the color data

    %% Load data
    disp("-->> Starting process");
    import templates.*;
    import guide.functions.split_hemisphere
    import guide.functions.tess_smooth
    import guide.functions.tess_hemisplit
    
    %% Load test data
    Cortex                  = load("templates/Cortex.mat");
    template                = load("templates/axes.mat");
    colorMap                = load("templates/mycolormap_brain_basic_conn.mat");
    currentAxes             = template.axes;
    hemis                   = []; % left or right or []
    [Cortex, iHideVert]     = split_hemisphere(Cortex, hemis);
    Source                  = J;  % Ensure J is defined in the calling script
    sources_iv = J;  % Input sources to be plotted
    sources_iv = sources_iv / max(sources_iv);  % Normalize the data
    smoothValue             = 0.2;
    SurfSmoothIterations    = 10;
    Vertices                = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);

    %% Plotting section
    patch(ax, ...
        'Faces', Cortex.Faces, ...
        'Vertices', Vertices, ...
        'FaceVertexCData', sources_iv, ...  % Use sources_iv for coloring
        'FaceColor', 'interp', ...
        'EdgeColor', 'none', ...
        'AlphaDataMapping', 'none', ...
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
    
    view(ax, 0, 0);  % Adjust view for 3D plots if necessary
    set(ax, 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w');
    
    % Directly apply the colormap and adjust color limits
    if exist('colorLimits', 'var') && ~isempty(colorLimits)
        caxis(ax, colorLimits);  % Set color axis limits
    else
        caxis(ax, 'auto');  % Automatically scale if no colorLimits are provided
    end
    
    % Apply the colormap here
    colormap(ax, 'parula');  % Can be changed depending on the requirement
    axis(ax, 'equal');
    rotate3d(ax, 'on');
    axis off;  % Turn off axes
end

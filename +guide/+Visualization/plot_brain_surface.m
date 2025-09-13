function plot_brain_views(J)
% plot_brain_views  Visualize cortical data (J) in three styled views
%   Layout: Large dorsal view (left), lateral & occipital stacked on right.
%
% INPUT:
%   J - data vector (Nv × 1), one value per cortical vertex

    disp("-->> Starting cortical surface visualization");

    %% === Imports and template loading ===
    import templates.*
    import guide.functions.split_hemisphere
    import guide.functions.tess_smooth
    import guide.functions.tess_hemisplit

    Cortex   = load("templates/Cortex.mat");
    colorMap = load("templates/mycolormap_brain_basic_conn.mat");

    % Hemisphere split
    hemis = [];
    [Cortex, iHideVert] = split_hemisphere(Cortex, hemis);

    %% === Prepare data ===
    J(iHideVert) = [];
    Vertices     = Cortex.Vertices;
    Vertices(iHideVert,:) = [];

    %% === Data smoothing (Gaussian) ===
    h = 0.004;   % bandwidth in mm
    D = pdist2(Vertices, Vertices);
    W = exp(-(D.^2) / (2*h^2));
    W = W ./ sum(W,2);

    J_smoothed = W * J;
    J_smoothed = J_smoothed * (norm(J(:)) / norm(J_smoothed(:)));

    %% === Geometry smoothing (mesh only) ===
    smoothValue      = 0.2;
    smoothIterations = 10;
    Vertices         = tess_smooth(Cortex.Vertices, ...
                                   smoothValue, ...
                                   smoothIterations, ...
                                   Cortex.VertConn, 1);

    %% === Custom layout with manual axes ===
    fig = figure('Color','w','Renderer','opengl','RendererMode','manual');

    % Axes positions: [x y width height] (normalized units 0–1)
    ax1 = axes('Position',[0.05 0.1 0.52 0.8]);  % dorsal view (slightly wider)
    ax2 = axes('Position',[0.60 0.55 0.35 0.35]); % lateral view (closer to left)
    ax3 = axes('Position',[0.60 0.10 0.35 0.35]); % occipital view


    % --- Dorsal ---
    plot_single_view(ax1, Cortex, Vertices, J_smoothed, colorMap.cmap_a, [0 90], 'Dorsal View');
    % --- Lateral ---
    plot_single_view(ax2, Cortex, Vertices, J_smoothed, colorMap.cmap_a, [90 0], 'Lateral View');
    % --- Occipital ---
    plot_single_view(ax3, Cortex, Vertices, J_smoothed, colorMap.cmap_a, [180 0], 'Occipital View');

    disp("-->> Visualization complete");
end

%% === Helper function ===
function plot_single_view(ax, Cortex, Vertices, data, cmap, viewAngles, labelText)
    patch(ax, ...
        'Faces',Cortex.Faces, ...
        'Vertices',Vertices, ...
        'FaceVertexCData',data, ...
        'FaceColor','interp', ...
        'EdgeColor','none', ...
        'FaceAlpha',0.99, ...
        'BackfaceLighting','lit', ...
        'AmbientStrength',0.2, ...
        'DiffuseStrength',0.5, ...
        'SpecularStrength',0.2, ...
        'SpecularExponent',1, ...
        'SpecularColorReflectance',0.5, ...
        'FaceLighting','gouraud', ...
        'EdgeLighting','gouraud');

    axis(ax,'equal','off');
    view(ax, viewAngles);
    colormap(ax, hot);
    title(ax, labelText, 'FontWeight','bold', 'FontSize',10, 'Color','k');
end

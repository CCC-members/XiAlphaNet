%%
%%  Description here
%%

%%
disp("-->> Starting process");


addpath("Visualization\ESI_plot\ESI_plot\app\");
addpath("Visualization\ESI_plot\ESI_plot\data\");
addpath("Visualization\ESI_plot\ESI_plot\functions\");
addpath("Visualization\ESI_plot\ESI_plot\template\");

%% Load test data
Cortex                  = load("tess_cortex_mid_high_8000V_fix");

template                = load("template/axes.mat");
colorMap                = load("template/mycolormap_brain_basic_conn.mat");
currentAxes             = template.axes;
hemis                   = []; % left or right or []
[Cortex, iHideVert]     = split_hemisphere(Cortex,hemis);
Source                  = J;%load("data/MEEG_source_alpha_7Hz_14Hz.mat");
J                       = J;%Source.J;
J(iHideVert)            = [];
sources_iv              = sqrt(abs(J));
sources_iv              = sources_iv/max(sources_iv(:));
smoothValue             = 0.66;
SurfSmoothIterations    = 10;
Vertices                = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);
fig = figure;
patch(currentAxes, ...
    'Faces',Cortex.Faces, ...
    'Vertices',Vertices, ...
    'FaceVertexCData',sources_iv, ...
    'FaceColor','interp', ...
    'EdgeColor','none', ...
    'AlphaDataMapping', 'none', ...
    'EdgeColor',        'none', ...
    'EdgeAlpha',        1, ...
    'BackfaceLighting', 'lit', ...
    'AmbientStrength',  0.5, ...
    'DiffuseStrength',  0.5, ...
    'SpecularStrength', 0.2, ...
    'SpecularExponent', 1, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting',     'gouraud', ...
    'EdgeLighting',     'gouraud', ...
    'FaceAlpha',.99);
set(currentAxes,'xcolor','w','ycolor','w','zcolor','w');
view(currentAxes,0,0);
colormap(currentAxes,colorMap.cmap_a);
set(currentAxes,"Parent",fig);
axis(currentAxes,'equal');
rotate3d(currentAxes,'on');







%%
disp("Process finished");
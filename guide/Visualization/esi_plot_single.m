%%
%%  Description here
%%

%%
disp("-->> Starting process");

%% Load test data
Cortex                  = load("templates/Cortex.mat");

template                = load("template/axes.mat");
colorMap                = load("template/mycolormap_brain_basic_conn.mat");
currentAxes             = template.axes;
hemis                   = []; % left or right or []
[Cortex, iHideVert]     = split_hemisphere(Cortex,hemis);
Source                  = J;

%% Original Code

% Original Code
J                       = J; %Source.J;
J(iHideVert)            = [];

% Apply smoothing to J using 'movmean' (moving average) or 'gaussian'
J_smoothed              = J;%smoothdata(J, 'movmean', 20);  % 5 is the window size, adjust as needed
%L = parameters.Model.L;
%sigma = 1;
%K = exp(-L.^2/(2*sigma^2));
%K = K./sum(K,2);
%J_smoothed = K*J_smoothed;
sources_iv = J_smoothed;%_star;
sources_iv              = sources_iv/max(sources_iv(:));
%%
smoothValue             = 0.4;
SurfSmoothIterations    = 10;
Vertices                = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);
fig = figure;
patch(currentAxes, ...
  'Faces',Cortex.Faces, ...
    'Vertices',Vertices, ...
    'FaceVertexCData',Cortex.SulciMap*0.00+1+sources_iv, ...
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
colormap("hot")
set(gcf,'Color','w');

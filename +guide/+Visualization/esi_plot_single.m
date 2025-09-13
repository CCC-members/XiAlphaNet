
disp("-->> Starting process");
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit
% Load test data
Cortex                  = load("templates/Cortex.mat");
template                = load("templates/axes.mat");
colorMap                = load("templates/mycolormap_brain_basic_conn.mat");
currentAxes             = template.axes;
hemis                   = []; % left or right or []
[Cortex, iHideVert]     = split_hemisphere(Cortex, hemis);
Source                  = J;
% Original Code
J                       = J; %Source.J;
J(iHideVert)            = [];

% --- Distance-based Gaussian smoothing ---
Vertices = Cortex.Vertices;
Vertices(iHideVert, :) = [];

% Parameters
h = 0.004; %0.004  % bandwidth (in mm), controls smoothing strength

% Compute Euclidean distance between vertices (NxN)
D = pdist2(Vertices, Vertices);  % warning: dense! for <5000 vertices

% Build Gaussian kernel
W = exp(-(D.^2) / (2 * h^2));
W = W ./ sum(W, 2);  % Normalize rows

% Smooth J
J_smoothed = W * J;
J_smoothed = J_smoothed * (norm(J(:)) / norm(J_smoothed(:)));


% Assign result
sources_iv = J_smoothed;

% Surface smoothing of vertices (geometry only)
smoothValue             = 0.5;
SurfSmoothIterations    = 20;
Vertices                = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);

% Plotting
fig = figure('Renderer','opengl');
set(gcf,'Renderer','zbuffer');           % stable 3D in software
set(gcf,'RendererMode','manual');        % lock the renderer
patch(currentAxes, ...
  'Faces',Cortex.Faces, ...
  'Vertices',Vertices, ...
  'FaceVertexCData',(sources_iv), ...
  'FaceColor','interp', ...
  'EdgeColor','none', ...
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
set(currentAxes,'xcolor','w','ycolor','w','zcolor','w');
view(currentAxes, 0, 0);
colormap(currentAxes, colorMap.cmap_a);
set(currentAxes, "Parent", fig);
axis(currentAxes, 'equal');
rotate3d(currentAxes, 'on');
%colormap("parula");
set(gcf, 'Color', 'w');



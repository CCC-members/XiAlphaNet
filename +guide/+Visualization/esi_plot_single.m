disp("-->> Starting process");

%% === Imports and template loading ===
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit

% Load template cortex & colormaps
Cortex   = load("templates/Cortex.mat");
template = load("templates/axes.mat");   % <-- template figure + axes
colorMap = load("templates/mycolormap_brain_basic_conn.mat");

% Get axes from template
fig         = figure('Renderer','opengl','Color','w');
currentAxes = copyobj(template.axes, fig);   % <-- copy axes into new figure
set(currentAxes,'Position',[0.05 0.05 0.9 0.9]); % adjust size

% Choose atlas
CortexiAtlas = 13;  % atlas to display borders
Atlas = Cortex.Atlas(CortexiAtlas);

% Split hemispheres
hemis = [];
[Cortex, iHideVert] = split_hemisphere(Cortex, hemis);

%% === Example data J (replace with yours) ===
J(isnan(J)) = 0;
J(iHideVert) = [];

%% === Smooth values ===
Vertices = Cortex.Vertices;
Vertices(iHideVert,:) = [];

h = 0.004;
D = pdist2(Vertices, Vertices);
W = exp(-(D.^2) / (2*h^2));
W = W ./ sum(W, 2);
J_smoothed = W * J;
J_smoothed = J_smoothed * (norm(J(:)) / norm(J_smoothed(:)));
sources_iv = J_smoothed;

% Smooth surface geometry
smoothValue          = 0.8;
SurfSmoothIterations = 20;
Vertices = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);

%% === Build atlas label per vertex ===
labels = zeros(size(Cortex.Vertices,1),1);
for r = 1:numel(Atlas.Scouts)
    labels(Atlas.Scouts(r).Vertices) = r;
end
labels(iHideVert) = 0;

%% === Plot cortical data on template axes ===
patch(currentAxes, ...
    'Faces',Cortex.Faces, ...
    'Vertices',Vertices, ...
    'FaceVertexCData',log(1+sources_iv), ...
    'FaceColor','interp', ...
    'EdgeColor','none', ...
    'BackfaceLighting','lit', ...
    'AmbientStrength',0.5, ...
    'DiffuseStrength',0.5, ...
    'SpecularStrength',0.2, ...
    'SpecularExponent',1, ...
    'SpecularColorReflectance',0.5, ...
    'FaceLighting','gouraud', ...
    'FaceAlpha',0.99);

axis(currentAxes,'equal','off');
view(currentAxes,[90 0]);   % default view
colormap(currentAxes,colorMap.cmap_a);
%colormap(pink)
rotate3d(fig,'on');

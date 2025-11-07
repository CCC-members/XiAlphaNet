%% === Directional Granger Difference (Permutation-based, GC already computed) ===

% === Load cortex and template ===
Cortex   = load("templates/Cortex_with_myelin.mat");
template = load("templates/axes.mat");
Cortex = Cortex.Cortex;
h =zscore(Cortex.Atlas(Cortex.iAtlas).MyelinValues);   % one hierarchy value per ROI

% === Input: precomputed Granger causality matrix ===
G = A;      % or Gxi_mean
N = size(G,1);
nPerm = 1000;
alpha_level = 0.002;

fprintf('Running %d permutations on GC for significance...\n', nPerm);

%% === Step 1: Permutation null on GC ===
% Null hypothesis: directed flow pattern is random w.r.t. ROI indices
G_null = zeros(N,N,nPerm);

for p = 1:nPerm
    perm = randperm(N);
    Gp = G(perm, perm);      % permute rows & columns jointly
    G_null(:,:,p) = Gp;
end

% --- Compute empirical p-values (one-tailed: large G)
pvals = zeros(N);
for i = 1:N
    for j = 1:N
        if i == j, continue; end
        null_dist = squeeze(G_null(i,j,:));
        pvals(i,j) = mean(null_dist >= G(i,j));
    end
end

mask_sig = pvals < alpha_level;
fprintf('Significant GC pairs: %.2f%% (%d of %d)\n', ...
    100*mean(mask_sig(:)), nnz(mask_sig), N^2);

%% === Step 2: Compute Directional Asymmetry Index (DAI) only for significant pairs ===
fprintf('Computing DAI for significant GC connections...\n');
eps_val = 1e-12;
DAI = zeros(N);

for i = 1:N
    for j = i+1:N
        %if mask_sig(i,j) && mask_sig(j,i)
            num = G(i,j) - G(j,i);
            den = G(i,j) + G(j,i) + eps_val;
            val = sign(h(i)-h(j))*num ;
            DAI(i,j) = val;
            DAI(j,i) = -val;
       % end
    end
end

% Keep only significant DAI values
DAI(~(mask_sig)) = 0;
if all(DAI(:) == 0)
    warning('No significant DAI connections found.');
end

% Normalize for visualization
DAI = DAI / max(abs(DAI(:)) + eps);
[src, tgt, val] = find(triu(DAI,1));
fprintf('Plotting %d significant directed DAI connections.\n', numel(src));



%scaledStrengths = scaledStrengths/max(abs(scaledStrengths));
% ROI coordinates
Scouts = Cortex.Atlas(13).Scouts;
connVertices = zeros(N,3);
for r = 1:N
    verts = Scouts(r).Vertices;
    connVertices(r,:) = Cortex.Vertices(verts(1),:);  % centroide 3D
end

%% === Step 5: Plot directional arrows ===
%% === Directional Asymmetry Plot (DAI ? [-1, 1]) ===
fprintf('Plotting %d DAI connections (-1 = feedback, +1 = feedforward)...\n', numel(src));
% === Create figure and copy template axes ===
fig = figure('Color', 'w', 'Renderer', 'opengl');
ax  = copyobj(template.axes, fig);
ax.Color = 'w';
hold(ax, 'on');

% === Plot cortex surface ===
patch(ax, 'Faces', Cortex.Faces, ...
          'Vertices', Cortex.Vertices, ...
          'FaceColor', [0.85 0.82 0.78], ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.25);
lighting(ax, 'gouraud');
camlight(ax, 'headlight');
view(ax, 3);

%% === HOT colormap with expanded orange/red midrange ===
nColors = 256;
baseMap = hot(nColors);

% Nonlinear remap: push more samples into midrange
% exponent <1 -> expands mid region; try 0.6-0.8
t = linspace(0,1,nColors)'.^2;

% Re-index base colormap with this nonlinear spacing
idx = round(1 + t*(nColors-1));
cmap = baseMap(idx,:);

% Warm tweak: cap yellow, keep strong orange
cmap(:,2) = cmap(:,2) * 0.85 + 0.10;   % soften yellow to ~0.8
cmap(:,3) = cmap(:,3) * 0.5;           % reduce blue -> deeper orange
cmap = max(min(cmap,1),0);

colormap(ax,cmap);
caxis(ax,[-1 1]);



%% === DAI scaling: compress extremes so 0.5<|DAI|<1 use more colors ===
beta = 2;    % lower beta = more compression; try 2–3
scaledStrengths = tanh(beta * val);   % replaces tanh(4*val)



% === Normalize and threshold strengths ===
scaledStrengths = max(min(scaledStrengths, 1), -1);
threshold = prctile(abs(scaledStrengths), 0);  % hide weakest 10%

% === Line width scaling ===
lineWidthMin = 0.3;
lineWidthMax = 6;

%% === Plot directional arrows ===
for k = 1:numel(src)
    val = scaledStrengths(k);
    if abs(val) < threshold, continue; end

    p1 = connVertices(src(k), :);
    p2 = connVertices(tgt(k), :);
    if any(~isfinite(p1)) || any(~isfinite(p2)), continue; end

    % Map [-1,1] ? [1,nColors]
    s = (val + 1) / 2;
    colorIdx = round(1 + s * (nColors - 1));
    colorIdx = max(min(colorIdx, nColors), 1);
    c = cmap(colorIdx, :);

    % Line width ? |value|
    lw = lineWidthMin + (lineWidthMax - lineWidthMin) * abs(val);

    dp = p2 - p1;
    quiver3(ax, p1(1), p1(2), p1(3), ...
                dp(1), dp(2), dp(3), ...
                0.9, ...
                'Color', c, ...
                'LineWidth', lw, ...
                'MaxHeadSize', 0.6);
end

% === Mark selected vertices ===
uniqueVertices = unique([src; tgt]);
for k = 1:numel(uniqueVertices)
    pos = connVertices(uniqueVertices(k), :);
    if all(isfinite(pos))
        plot3(ax, pos(1), pos(2), pos(3), 'o', ...
              'MarkerSize', 6, ...
              'MarkerEdgeColor', 'k', ...
              'MarkerFaceColor', 'w', ...
              'LineWidth', 1.5);
    end
end

axis(ax, 'equal', 'vis3d', 'off');

% === Title and Colorbar ===
title(ax, ...
    'Significant GC  Directional Asymmetry Index (-1 = Feedback, +1 = Feedforward)', ...
    'FontWeight', 'bold', ...
    'FontSize', 13);

cb = colorbar(ax, 'Location', 'eastoutside');
cb.Label.String = sprintf('DAI (\\gamma = %.2f)', beta);
cb.Limits = [-1 1];
cb.Ticks = -1:0.5:1;
cb.TickLabels = {'-1','-0.5','0','0.5','1'};

fprintf('--- DAI cortical plot complete ---\n');

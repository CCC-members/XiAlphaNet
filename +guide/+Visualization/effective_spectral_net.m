%% === Load cortex surface and ROI information ===
Cortex = load("templates/Cortex.mat");
template = load("templates/axes.mat");

% === Input: Cross-spectrum matrix A (360x360 complex) ===
N = 360;

% === Step 1: Threshold based on raw cross-spectrum magnitude ===
rawStrength = abs(A);
threshold = prctile(rawStrength(:), 99.0);  % Top 1%
mask = triu(rawStrength > threshold, 1);
[src, tgt] = find(mask);
fprintf('Selected %d connections from raw |A| > %.4f\n', numel(src), threshold);

if isempty(src)
    warning('No connections meet raw cross-spectrum threshold!');
    return
end

% === Step 2: Compute isolated partial coherence ===
% Input: A is the cross-spectrum matrix (complex, NxN)
%        src, tgt are indices of connections (i.e., i <- j)
P = pinv(A);  % Precision matrix (complex)
Granger = zeros(N);  % Granger causality matrix

for k = 1:numel(src)
    i = src(k); j = tgt(k);
    
    % Elements of the precision matrix
    pii = real(P(i,i));
    pjj = real(P(j,j));
    pij = real(P(i,j));

    % Conditional spectral power of i given j is removed
    S_ii_cond = 1 / (pii - (pij^2 / pjj));
    % Full spectral power of i
    S_ii_full = 1 / pii;

    % Granger causality i <- j
    if S_ii_cond > 0 && S_ii_full > 0
        Granger(i,j) = log(S_ii_full / S_ii_cond);
    end
end


% === Step 3: Normalize and apply power-law contrast mapping ===
strengths = PCoh(sub2ind(size(Granger), src, tgt));
normStrengths = strengths / max(strengths);  % Normalize to [0, 1]

gamma = 0.3;  % Power-law contrast
scaledStrengths = normStrengths .^ gamma;

% === Get ROI centroid vertex ===
Scouts = Cortex.Atlas(8).Scouts;
connVertices = zeros(360, 3);
for r = 1:360
    verts = Scouts(r).Vertices;
    connVertices(r, :) = Cortex.Vertices(verts(1), :);  % First vertex
end

% === Create figure and copy template axes ===
fig = figure('Color', 'w', 'Renderer', 'opengl');
ax = copyobj(template.axes, fig);
ax.Color = 'w';
hold(ax, 'on');

% === Plot cortex surface ===
patch(ax, 'Faces', Cortex.Faces, ...
          'Vertices', Cortex.Vertices, ...
          'FaceColor', [0.8 0.78 0.75], ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.2);
lighting(ax, 'gouraud');
camlight(ax, 'headlight');
view(ax, 3);

% === Step 4: Custom Hot-to-YellowOrange Colormap ===
nColors = 256;
x = linspace(0, 1, nColors)';
r = min(1, 2 * x + 0.2);
g = min(1, 1.5 * max(x - 0.2, 0));
b = min(1, 1.2 * max(x - 0.7, 0)) * 0.2;
cmap = [r, g, b];
backgroundColor = [0.1 0 0];  % Deep red/black background base

% === Step 5: Plot directional arrows ===
threshold = prctile(scaledStrengths, 10);  % Hide weakest 10%

for k = 1:numel(src)
    if scaledStrengths(k) < threshold
        continue
    end

    p1 = connVertices(src(k), :);
    p2 = connVertices(tgt(k), :);
    if any(~isfinite(p1)) || any(~isfinite(p2)), continue; end

    s = (scaledStrengths(k) - threshold) / (1 - threshold);
    s = max(min(s, 1), 0);
    colorIdx = round(1 + 255 * s);
    colorIdx = max(min(colorIdx, 256), 1);

    alphaVal = 0.3 + 0.7 * s;
    baseColor = cmap(colorIdx, :);
    blendedColor = alphaVal * baseColor + (1 - alphaVal) * backgroundColor;

    % Arrow from source to target
    dp = p2 - p1;
    quiver3(ax, p1(1), p1(2), p1(3), ...
                dp(1), dp(2), dp(3), ...
                0.9, ...
                'Color', blendedColor, ...
                'LineWidth', 2 + 4 * s, ...   % Wider arrows
                'MaxHeadSize', 0.6);
end

% === Step 6: Mark selected vertices ===
uniqueVertices = unique([src; tgt]);
for k = 1:numel(uniqueVertices)
    idx = uniqueVertices(k);
    pos = connVertices(idx, :);
    if all(isfinite(pos))
        plot3(ax, pos(1), pos(2), pos(3), 'o', ...
              'MarkerSize', 6, ...
              'MarkerEdgeColor', 'k', ...
              'MarkerFaceColor', 'w', ...
              'LineWidth', 1.5);
    end
end

axis(ax, 'equal', 'vis3d', 'off');
title(ax, sprintf('Isolated Partial Coherence (Power Contrast, ? = %.2f)', gamma));

% === Step 7: Add colorbar ===
colormap(ax, cmap);
cb = colorbar(ax, 'Location', 'eastoutside');
cb.Label.String = sprintf('Partial Coherence (Power Transform ? = %.2f)', gamma);
cb.Ticks = [0, 0.25, 0.5, 0.75, 1];
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x.^(1/gamma)), cb.Ticks, 'UniformOutput', false);

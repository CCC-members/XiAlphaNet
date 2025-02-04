function circularGraphWithWeights(CM, varargin)
    % CM is the connectivity matrix with log-scaled weights

    % Check the number of inputs
    p = inputParser;
    addParameter(p, 'Colormap', 'hot', @isstring);
    parse(p, varargin{:});
    cmap = colormap(p.Results.Colormap);  % Use the specified colormap

    % Number of nodes
    N = size(CM, 1);

    % Create figure and axis
    figure('Position', [300, 300, 560, 560], 'Color', [1, 1, 1]);
    axes('Position', [0.1, 0.1, 0.8, 0.8], 'Color', [1, 1, 1]);

    % Angle for each node
    theta = linspace(0, 2*pi, N+1)';
    theta(end) = [];

    % Node positions
    x = cos(theta);
    y = sin(theta);

    % Draw the edges
    hold on;
    for i = 1:N
        for j = 1:N
            if CM(i, j) > 0
                % Weight to color mapping
                w = CM(i, j);
                colorIdx = min(round(w * size(cmap, 1)), size(cmap, 1));  % Ensure index is within colormap bounds
                thisColor = cmap(colorIdx, :);

                % Draw edge
                plot([x(i) x(j)], [y(i) y(j)], 'LineWidth', 2, 'Color', thisColor);
            end
        end
    end

    % Draw the nodes
    for i = 1:N
        plot(x(i), y(i), 'o', 'Color', [0, 0, 0], 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerSize', 10);
    end

    hold off;
    axis equal off;
end

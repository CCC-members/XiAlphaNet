function threshold = set_threshold_em(V)
    % Remove NaN values and reshape the data
    V = V(~isnan(V));
    V = V(:);

    % Filter out near-zero values
    V = V(abs(V) > 1e-10);

    % Fit a two-component GMM using the EM algorithm
    try
        gm = fitgmdist(V, 2, 'Replicates', 10, ...
            'Options', statset('MaxIter', 500, 'TolFun', 1e-6));
    catch ME
        % If GMM fitting fails, use the 90th percentile as a fallback
        threshold = prctile(V, 90);
        return;
    end

    % Extract parameters for the two components
    mu = gm.mu;
    sigma = gm.Sigma;
    weights = gm.ComponentProportion;

    % Sort the components by their means
    [~, idx] = sort(mu);
    mu_sorted = mu(idx);
    sigma_sorted = sigma(:, :, idx);
    weights_sorted = weights(idx);

    % For 1D data, extract the variances (diagonal elements)
    sigma_1 = sigma_sorted(1, 1, 1);
    sigma_2 = sigma_sorted(1, 1, 2);

    % Define the probability density functions for each Gaussian component
    pdf_1 = @(x) (1 / sqrt(2 * pi * sigma_1)) * exp(-(x - mu_sorted(1)).^2 / (2 * sigma_1));
    pdf_2 = @(x) (1 / sqrt(2 * pi * sigma_2)) * exp(-(x - mu_sorted(2)).^2 / (2 * sigma_2));

    % Define the intersection function
    intersection_func = @(x) pdf_1(x) - pdf_2(x);

    % Start with an interval multiplier k = 0
    k = 0;
    % Initial endpoints: starting from the means (k=0)
    a = mu_sorted(1) - k * sigma_1;
    b = mu_sorted(2) + k * sigma_2;

    % Set maximum iterations for expanding the interval
    maxIter = 20;
    iter = 0;
    
    % Expand the interval until the endpoints bracket a root (i.e. f(a)*f(b)<0)
    while iter < maxIter && intersection_func(a)*intersection_func(b) >= 0
        k = k + 1;
        a = mu_sorted(1) - k * sigma_1;
        b = mu_sorted(2) + k * sigma_2;
        iter = iter + 1;
    end

    if intersection_func(a)*intersection_func(b) < 0
        % Found a valid interval, so use fzero
        threshold = fzero(intersection_func, [a, b]);
    else
        % If a proper bracket wasn't found after maxIter, fallback to the 90th percentile
        threshold = prctile(V, 90);
    end

    % Optional: Plotting for visualization (remove in production)
    %{
    test_vals = linspace(min(V), max(V), 1000);
    figure;
    hold on;
    plot(test_vals, pdf_1(test_vals), 'r-', 'LineWidth', 1.5);
    plot(test_vals, pdf_2(test_vals), 'g-', 'LineWidth', 1.5);
    plot(test_vals, intersection_func(test_vals), 'b-', 'LineWidth', 1.5);
    plot(threshold, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Value');
    ylabel('Probability Density / Difference');
    title('Intersection of Gaussian Components');
    grid on;
    hold off;
    %}
end

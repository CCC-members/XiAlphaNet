function threshold = set_threshold_em(V)
    % This function sets a threshold for deleting near-zero values of a sparse vector
    % using an Expectation-Maximization (EM) based Gaussian Mixture Model (GMM).
    % The threshold is where the two Gaussian distributions intersect.

    % Step 1: Remove NaN values and reshape the data
    V = V(~isnan(V));  % Remove NaN values
    V = V(:);          % Ensure it's a column vector

    % Step 2: Remove very small values or zeros if necessary to improve fitting
    V = V(abs(V) > 1e-6);  % Filter out near-zero values (adjust threshold if needed)


    % Step 4: Fit a Gaussian Mixture Model with two components (EM algorithm)
    try
        gm = fitgmdist(V, 2, 'Replicates', 10, 'Options', statset('MaxIter', 500, 'TolFun', 1e-6));
    catch ME
        %warning('GMM fitting failed: %s', ME.message);
        threshold = prctile(V,90);  % Return NaN if fitting fails
        return;
    end

    % Step 5: Get the parameters of the two Gaussian distributions
    mu = gm.mu;                    % Means of the two components
    sigma = gm.Sigma;              % Covariances of the two components
    weights = gm.ComponentProportion; % Component weights

    % Step 6: Identify the Gaussian components corresponding to near-zero and non-zero
    [~, idx] = sort(mu);           % Sort the means to identify the zero and non-zero components
    mu_sorted = mu(idx);
    sigma_sorted = sigma(:, :, idx);  % Sort the covariances
    weights_sorted = weights(idx);

    % Step 7: For 1D data, extract the variances (diagonal elements) from sigma_sorted
    sigma_1 = sigma_sorted(1, 1, 1);  % Variance of the first Gaussian
    sigma_2 = sigma_sorted(1, 1, 2);  % Variance of the second Gaussian

    % Step 8: Define the probability density functions for the two Gaussian components
    pdf_1 = @(x) (1 / sqrt(2 * pi * sigma_1)) * exp(-(x - mu_sorted(1)).^2 / (2 * sigma_1));
    pdf_2 = @(x)  (1 / sqrt(2 * pi * sigma_2)) * exp(-(x - mu_sorted(2)).^2 / (2 * sigma_2));

    % Step 9: Define the intersection function
   test_vals = linspace(min(min(V),0),max(V),1000);
    intersection_func = @(x) pdf_1(x) - pdf_2(x);
    f_vals = intersection_func (test_vals);
 

    % Step 12: Find the threshold where the intersection_func is closest to zero
    threshold = fzero(intersection_func,[mu_sorted(1),mu_sorted(2)]);

    % % %%Optional: Plotting for visualization (can be removed in production)
    % figure;
    % hold on 
    % p1 = pdf_1(test_vals);
    % p2 = pdf_2(test_vals);
    % plot(test_vals,p1,'r-', 'LineWidth', 1.5)
    % plot(test_vals,p2,'g-', 'LineWidth', 1.5)
    % plot(test_vals, f_vals, 'b-', 'LineWidth', 1.5);
    % hold on;
    % plot(threshold, f_vals(idx_min), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    % xlabel('Value');
    % ylabel('PDF Difference');
    % title('Intersection of Gaussian Components');
    % grid on;
    % hold off;
end

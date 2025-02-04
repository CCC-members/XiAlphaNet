% Function to generate Gaussian activation based on distance from selected vertices
function activation = generateGaussianActivation(selected_vertices, L, mean_value, sigma_spread, Nr)
    % Preallocate the activation vector
    activation = zeros(Nr, 1);
    
    % Loop over each vertex and compute the activation based on the distance matrix L
    for i = 1:Nr
        % Calculate distance from selected vertices to the i-th vertex
        dists = L(i, selected_vertices);  % Distances between i-th vertex and all selected vertices
        min_dist = min(dists);  % Choose the minimum distance (closest point in the region)
        
        % Generate Gaussian activation based on distance
        activation(i) = mean_value * exp(-min_dist^2 / (2 * sigma_spread^2));
    end
end

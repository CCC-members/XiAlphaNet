% Function to generate Gaussian activation symmetrically across hemispheres for a(:,1)
function activation = generateGaussianActivationSymmetric(verticesL, verticesR, L, mean_value, sigma_spread, Nr)
    % Preallocate the activation vector
    activation = zeros(Nr, 1);
    
    % Generate Gaussian activation for the left hemisphere
    activationL = generateGaussianActivation(verticesL, L, mean_value, sigma_spread, Nr);
    
    % Symmetrize activation for the right hemisphere
    activationR = generateGaussianActivation(verticesR, L, mean_value, sigma_spread, Nr);
    
    % Combine left and right activations symmetrically
    activation = (activationL + activationR) / 2;
end

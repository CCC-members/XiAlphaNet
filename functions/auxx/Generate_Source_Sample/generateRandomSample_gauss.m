function x = generateRandomSample_gauss(parameters, var)
    Nr = parameters.Dimensions.Nr;
    t = 0.0067;%activation_threshold(parameters);
    L = parameters.Model.L;
    
    % Set the target means for activations
    meanE1 = (5/6) * t;
    meanE2 = 0.001;
    meanE3 = 2.5;
    
    meanA1 = (1/6) * t;
    meanA2 = 0.008;
    meanA3 = 15;
    meanA4 = 10;  % Midpoint of the range [5, 13]
    
    % Generate e(:,1) quasi-uniform across the brain and symmetric
    verticesL = getRandomVerticesInRegion('left');  % Get random left hemisphere vertices
    verticesR = getRandomVerticesInRegion('right'); % Get random right hemisphere vertices
    spread_e1 = 100;%20+80*rand(1);
    spread_a1 = 10*rand(1);
    e1 = generateGaussianActivationSymmetric(verticesL, verticesR, L, meanE1, spread_e1, Nr);
    
    % Generate a(:,1) restricted to occipital and parietal areas, symmetric
    verticesL = getRandomVerticesInRegion('left');  % Get random left hemisphere vertices
    verticesR = getRandomVerticesInRegion('right'); % Get random right hemisphere vertices
    a1 = generateGaussianActivationSymmetric(verticesL, verticesR, L, meanA1, spread_a1, Nr);  % Symmetric Gaussian
    
    % Generate the remaining columns for e and a
    e2 = abs(normrnd(meanE2, sqrt(var)/1000, Nr, 1));
    e3 = gamrnd(meanE3^2 / var, var / meanE3, Nr, 1);
    e = [e1, e2, e3];
    
    a2 = abs(normrnd(meanA2, sqrt(var)/1000, Nr, 1));
    a3 = gamrnd(meanA3^2 / var, var / meanA3, Nr, 1);
    a4 = normrnd(meanA4, sqrt(var), Nr, 1);
    a4 = max(min(a4, 13), 5);  % Clip a4 to the range [5, 13]
    a = [a1, a2, a3, a4];
    
    % Reshape matrices e and a into vectors and concatenate with sigma^2
    sigma2 = gamrnd(1, 1);  % Shape=10, Rate=1, Mean=10
    x = v2x(e, a, sigma2);
end

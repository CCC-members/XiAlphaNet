% Parameters for the Gamma distribution
shape = 2;  % Example shape parameter
scale = 2;  % Example scale parameter

% Generate the initial vector of 8000 Gamma-distributed random numbers
vec = gamrnd(shape, scale, [8003, 1]);

% Define the indices for the segments that need to be larger or equal
indices_larger = [1:1900, 3700:4900, 7500:8000];

% Increase values at specific indices by adding a constant or scaling
increment = max(vec) * 0.2;  % Calculate an increment as 20% of the max value in the vector
vec(indices_larger) = vec(indices_larger) + increment;

% Optionally, ensure that these values remain Gamma distributed by scaling:
% This is a rough approximation and may not perfectly align with a true Gamma distribution
scale_factor = 1.2;  % Scaling factor to adjust the distribution slightly
vec(indices_larger) = vec(indices_larger) * scale_factor;

% Display or use the vector
disp(vec);
J = vec;
esi_plot;
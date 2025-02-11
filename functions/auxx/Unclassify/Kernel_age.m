function K = Kernel_age(age_star, age, h)
    % Kernel_age - Compute the Gaussian kernel for age differences.
    %
    % Inputs:
    %   age_star - The target age (scalar or vector).
    %   age - A vector of ages to compare with (can be the same size as age_star).
    %   h - The bandwidth parameter for the Gaussian kernel.
    %
    % Output:
    %   K - The Gaussian kernel values (same size as age).

    % Compute the difference between age_star and age
    age_diff = age_star - age;
    
    % Gaussian kernel formula
    K = exp(-(age_diff).^2 / (2 * h^2));
end

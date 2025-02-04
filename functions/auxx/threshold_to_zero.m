function A = threshold_to_zero(A)
    % Find the maximum value in the vector A
    maxVal = max(A);
    
    % Calculate the threshold (10% of the maximum value)
    threshold = 0.1 * maxVal;
    
    % Set values less than the threshold to zero
    A(A < threshold) = 0;
end
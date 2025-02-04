function M = v2m(v)
    % v2m - Inverse half-vectorization to form a symmetric matrix.
    %
    % Inputs:
    %   v - A vector containing the upper (or lower) triangular part of M.
    %
    % Output:
    %   M - The symmetric matrix (n x n) reconstructed from v.

    % Calculate the size of the original matrix
    n = (-1 + sqrt(1 + 8 * length(v))) / 2;
    
    % Ensure the length of v corresponds to a valid symmetric matrix
    if mod(n, 1) ~= 0
        error('Length of vector does not correspond to a valid symmetric matrix.');
    end
    
    % Initialize the matrix
    M = zeros(n);

    % Fill the upper triangular part (including the diagonal) with the vector elements
    M(triu(true(n))) = v;
    
    % Mirror the upper triangular part to the lower triangular part
    M = M + triu(M, 1)';
end

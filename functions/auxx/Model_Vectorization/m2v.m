function v = m2v(M)
    % m2v - Half-vectorize a symmetric matrix M.
    %
    % Inputs:
    %   M - A symmetric matrix (n x n).
    %
    % Output:
    %   v - A vector containing the upper (or lower) triangular part of M.

    % Ensure M is symmetric
    if ~issymmetric(M)
        error('Input matrix must be symmetric.');
    end
    
    % Get the upper triangular part of the matrix, including the diagonal
    v = M(triu(true(size(M))));
end

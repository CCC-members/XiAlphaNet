function v = m2v(M)
    % m2v - Half-vectorize a symmetric matrix M.
    %
    % Inputs:
    %   M - A symmetric matrix (n x n).
    %
    % Output:
    %   v - A vector containing the upper (or lower) triangular part of M.
    M = (M+M')/2;
    % Get the upper triangular part of the matrix, including the diagonal
    v = M(triu(true(size(M))));
end

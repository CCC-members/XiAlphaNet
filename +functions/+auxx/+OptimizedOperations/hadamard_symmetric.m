function C = hadamard_symmetric(A, B)
    % Validate that A and B are symmetric
    A = sym_matrix(A);
    B = sym_matrix(B);
    
    % Get the size of the matrices
    n = size(A, 1);

    % Initialize the result matrix C
    C = zeros(n);

    % Compute the Hadamard product only for the upper triangular part
    for i = 1:n
        for j = i:n  % Loop over the upper triangular part, including diagonal
            C(i, j) = A(i, j) * B(i, j);
            if i ~= j
                C(j, i) = C(i, j);  % Reflect the value to the lower triangular part
            end
        end
    end
end


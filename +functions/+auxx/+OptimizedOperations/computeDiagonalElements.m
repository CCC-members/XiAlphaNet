function diag_elements = computeDiagonalElements(T, D)
    n = size(T, 1);
    diag_elements = zeros(n, 1);
    TD = T * D;
    for i = 1:n
        diag_elements(i) = TD(i, :) * T(i, :)';
    end
end

function TDT = computeTDT(T, D)
    % Compute TDT' where D is the diagonal elements of a diagonal matrix.
    % T is an n-by-m matrix.
    % D is a vector containing the diagonal elements of the diagonal matrix.

    % Scale each column of T by the corresponding diagonal entry in D
    TD = T .* D';  % Element-wise multiplication to scale columns

    % Compute the product of the scaled matrix TD with T'
    TDT = TD * T';
    TDT = 0.5*(TDT+TDT');
end

function [A] = nulldiag(A);
[n,~] = size(A);
for j=1:n
    A(j,j) = 0;
end
end

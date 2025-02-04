function [g3] = evalg3(x,parameters)
% This function eval the regularition term of the objective function 

%Parameters
lambda3 = parameters.Regularization(3);


 % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[~,a,~]  =x2v(x);

% Calculate the number of rows in a
numRowsA = size(a, 1);

% Compute the second term involving norms of rows of e and a
g3 = lambda3 * sum(sum(sqrt(sum(a.^2, 2))))/numRowsA;

end

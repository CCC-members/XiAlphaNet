function [g2] = evalg2(x,parameters)
% This function eval the regularition term of the objective function 

%Parameters
lambda2 = parameters.Regularization(2);


 % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[e,~,~]  =x2v(x);

% Calculate the number of rows in e
numRowsE = size(e, 1);

% Compute the second term involving norms of rows of e and a
g2 = lambda2 * sum(sum(sqrt(sum(e.^2, 2))))/numRowsE;

end

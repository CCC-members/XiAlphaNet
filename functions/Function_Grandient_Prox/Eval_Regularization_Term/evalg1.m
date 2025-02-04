function [g1] = evalg1(x,parameters)
% This function eval the regularition term of the objective function 

%Parameters
lambda1 = parameters.Regularization(1);


 % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[~,~,sigma2]  =x2v(x);

% Compute the first term involving sigma^2
g1 = lambda1 * norm(sigma2)^2;

end

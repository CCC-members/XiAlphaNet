function [g1] = evalg1(x)
% This function eval the regularition term of the objective function 

 % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[~,~,sigma2]  =x2v(x);

% Compute the first term involving sigma^2
g1 = norm(sigma2)^2;

end

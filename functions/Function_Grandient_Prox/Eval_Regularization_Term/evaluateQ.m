function [term1,term2, qval] = evaluateQ(x,parameters)
% This function eval the regularition term of the objective function 

%Parameters
lambda1 = parameters.Regularization(1);
lambda2 = parameters.Regularization(2);

 % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[e,a,sigma2]  =x2v(x);

% Compute the first term involving sigma^2
term1 = lambda1 * norm(sigma2)^2;

% Compute the second term involving norms of rows of e and a
term2 = lambda2 * sum(sum(sqrt(sum(e.^2, 2))) + sum(sqrt(sum(a.^2, 2))));

% Evaluation
qval  = term1 + term2;
end

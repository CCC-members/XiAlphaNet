function [g] = evalg(x,opt)
import functions.*
import functions.auxx.*
import functions.auxx.ModelVectorization.*
% This function eval the regularition term of the objective function
switch opt
    case 1
        % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
        [~,~,sigma2]  =x2v(x);
        % Compute the first term involving sigma^2
        g = norm(sigma2)^2;
    case 2
        % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
        [e,~,~]  =x2v(x);
        % Calculate the number of rows in e
        numRowsE = size(e, 1);
        % Compute the second term involving norms of rows of e and a
        g = sum(sum(sqrt(sum(e.^2, 2))))/numRowsE;
    case 3
        % Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
        [~,a,~]  =x2v(x);
        % Calculate the number of rows in a
        numRowsA = size(a, 1);
        % Compute the second term involving norms of rows of e and a
        g = sum(sum(sqrt(sum(a.^2, 2))))/numRowsA;
end
end

function [n] = tensor_norm(X,p)
% Compute the normalized tensor norm of X with oreder p
dim = size(X);
N = prod(dim);
X = reshape(X,[1,N]);
n = norm(X,p)/N;
end

    
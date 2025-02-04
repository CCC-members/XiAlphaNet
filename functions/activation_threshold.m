function [t] = activation_threshold(parameters)
S = parameters.Data.Cross;
T = parameters.Model.T;
[~,~,f] =size(S);
tn = 0;
td = 0;
for j  = 1:f
    tn = tn + trace(S(:,:,j));
    td = td + trace(T(:,:,j)*T(:,:,j)');
end
t = tn/(td);
end

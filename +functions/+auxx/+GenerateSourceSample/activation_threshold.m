function [t] = activation_threshold(S,T)
[~,~,f] =size(S);
tn = 0;
td = 0;
for j  = 1:f
    tn = tn + real(trace(S(:,:,j)));
    td = td + real(trace(T(:,:,j)*T(:,:,j)'));
end
t = tn/(td);
end
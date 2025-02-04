function [scaled_x] = global_scale_factor(x,parameters)
S = parameters.Data.Cross;
[Nc,~,Nw]  = size(S);

% Initialize the variable for k
k = 0;
% Compute k based on the provided formula
for j = 1:Nw
    k = k + real(trace(diag(log(diag(S(:,:,j))))));
end
% Compute Gs
k = k / (Nc * Nw);
Gs = exp(k);
[e,a,s2] = x2v(x.Solution);
e(:,1) = e(:,1)/Gs;
a(:,1) = a(:,1)/Gs;
s2 = s2/Gs;
scaled_x = x;
scaled_x.Solution = v2x(e,a,s2);
end

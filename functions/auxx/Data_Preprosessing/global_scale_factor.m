function [S] = global_scale_factor(Cross)
S = Cross;
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
S = S/Gs;
end

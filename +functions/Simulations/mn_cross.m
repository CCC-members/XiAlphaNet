function [Sjj]=mn_cross(Svv,parameters)
K=parameters.Model.K;
K =pinv(K);
U_map =K;
Nw =parameters.Dimensions.Nw;
for i=1:Nw
    Sjj(:,:,i) = U_map*Svv(:,:,i)*U_map';
end
end

function [Sjj]=mn_cross(Svv,K)
[~,~,Nw] = size(Svv);
[~,~,Nk] = size(K);
%K = reshape(K,[Ne,Nw,Nk]);
if Nk == 1
    U_map = pinv(K);
else
    for i = 1:Nw
        U_map(:,:,i) = pinv(K(:,:,i));
    end
end

parfor i=1:Nw
    if Nk ==1
        Sjj(:,:,i) = U_map*Svv(:,:,i)*U_map';
    else
        Sjj(:,:,i) = U_map(:,:,i)*Svv(:,:,i)*U_map(:,:,i)';
    end
end
end

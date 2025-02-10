function [Svv] = remove_hemisphere_effect(Svv,Lvj,ind)
p             = size(Lvj,1);
scaleLvj      = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj           = Lvj/scaleLvj;
scaleV        = (sum(abs(diag(Svv)))/p);
Svv           = Svv/scaleV;
[U,D,V]       = svd(Lvj,'econ');
d             = diag(D);
gamma         = 1;
gammad2       = gamma*d.^2;
dvv           = gammad2./(gammad2 + 1);
Tjv           = V*diag(dvv)*U';
LvjTjv_ind    = Lvj(:,ind)*Tjv(ind,:);
Svv           = Svv - LvjTjv_ind*Svv*LvjTjv_ind';
Svv           = (Svv + Svv')/2;
Svv           = Svv*scaleV;
end
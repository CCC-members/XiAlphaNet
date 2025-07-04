function source=inverse(Svv,Lvj)

% input Lvj:  NxMxP leadfield tensor for N channels, M voxels, and 
%           P dipole directions. Typically P=3. (If you do MEG for 
%           a spherical volume conductor or reduce the rank, you must 
%           reduce L such that it has full rank for each voxel, such that,
%           e.g., P=2)


% [A,Wout]=mkfilt_eloreta(L,regu,W)

[source.eloreata.Sjj,source.eloreata.Tjv]=gcv_eloreta(Svv,Lvj);
[source.lcmv.Sjj,source.lcmv.Tjv]=gcv_lcmv(Svv,Lvj);

end



function [Sjj,Tjv]=gcv_eloreta(Svv,Lvj)
    p             = size(Lvj,1);
    Ip            = eye(p);
    param.gamma1          = 0.01;
    param.gamma2          = 0.5;
    param.delta_gamma     = 0.01;
    gamma1        = param.gamma1;
    gamma2        = param.gamma2;
    delta_gamma   = param.delta_gamma;
    gamma_grid    = gamma1:delta_gamma:gamma2;
    count         = 1;
    for gamma = gamma_grid
        [Tjv,Wout]       = mkfilt_eloreta(Lvj,gamma);
         Tjv              = Tjv';
         Txiv             = Ip - Lvj*Tjv;
         gcv(count)       = (1/p)*sum(abs(diag(Txiv*Svv*Txiv')))/((1/p)*sum(abs(diag(Txiv))))^2;
         count             = count + 1;
     end
     [gcv_opt,idx_gamma]       = min(gcv);
    gamma                     = gamma_grid(idx_gamma);
    [Tjv,Wout]                = mkfilt_eloreta(Lvj,gamma);
    Sjj=Tjv'*Svv*Tjv;
end



function [Sjj,Tjv]=gcv_lcmv(Svv,Lvj)
    p             = size(Lvj,1);
    Ip            = eye(p);
    param.gamma   = sum(abs(diag(Svv)))/(length(Svv)*100);
    p             = size(Lvj,1);
    scaleLvj      = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
    Lvj           = Lvj/scaleLvj;
    scaleV        = (sum(abs(diag(Svv)))/p);
    Svv           = Svv/scaleV;
    gamma         = param.gamma;
    gamma         = gamma/scaleV;
    [Tjv]   = mkfilt_lcmv(Lvj,Svv,gamma);
    Sjj=Tjv'*Svv*Tjv;
end









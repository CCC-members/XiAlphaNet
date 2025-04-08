function [measures] = quality_measures(sol,Theta0_sim)
Nsim     = size(sol,2);
Nmeth    = size(sol{1},3);
measures = zeros(Nmeth,5,Nsim);
for sim = 1:Nsim
%     Xsim = Theta_sim{sim}(:,:,numSlice);
    Xsim = Theta0_sim{sim};
%     Xsim = abs(Xsim)./(sqrt(abs(diag(Xsim))*abs(diag(Xsim))') + 1E-4*max(abs((diag(Xsim)))));
    Xsim = Xsim - diag(diag(Xsim));
    Xsim(abs(Xsim) > 0) = 1;
    Xest = sol{3,sim};
    for meth = 1:Nmeth
        Xtmp = Xest(:,:,meth);
%         Xtmp = abs(Xtmp)./(sqrt(abs(diag(Xtmp))*abs(diag(Xtmp))') + 1E-4*max(abs((diag(Xtmp)))));
        Xtmp = Xtmp - diag(diag(Xtmp));
        if max(abs(Xtmp(:))) ~= 0 
            Xtmp = abs(Xtmp)/max(abs(Xtmp(:)));
        end
        %% ROC analysis
        %X->spec Y->sens
        [spec,sens,th,auc,OPTROCPT] = perfcurve(Xsim(:),Xtmp(:),true);
        optINDX = find((spec == OPTROCPT(1))&(sens == OPTROCPT(2)));
        ppv = perfcurve(abs(Xsim(:)),abs(Xtmp(:)),true,'Xcrit','ppv');
        measures(meth,1,sim) = auc;
        measures(meth,2,sim) = sens(optINDX);
        measures(meth,3,sim) = 1 - spec(optINDX);
        measures(meth,4,sim) = ppv(optINDX);
        measures(meth,5,sim) = 2*ppv(optINDX)*sens(optINDX)/(ppv(optINDX)+sens(optINDX));
    end
end
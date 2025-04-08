function Sjj_est=call_higgs(Svv,Lvj)
Nsegments = 100;
penalty             = [1 2 0]; % 1 (lasso) 2 (frobenious) 0 (naive)
param.maxiter_outer = 60;
param.maxiter_inner = 30;
p                   = size(Svv,1);%sensor number
q                   = size(Lvj,2);
param.p             = p;
param.q             = q;
param.Ip            = eye(p);
param.Iq            = eye(q);
param.m             = Nsegments;
aj                  = sqrt(log(q)/Nsegments);                                           
Ajj_diag            = 0;                                                      
Ajj_ndiag           = 1;                                               
Ajj                 = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj            = aj;
param.Ajj           = Ajj;
param.axi           = 1E-4;
param.Axixi         = eye(p);
param.Axixi_inv     = eye(p);
param.ntry          = 0;
param.prew          = 0;
param.nu            = Nsegments;
param.rth1          = 0.7;
param.rth2          = 3.16;
run_bash_mode = 'true'
param.run_bash_mode   = run_bash_mode;
param.Op              = ones(p,1);
param.run_bash_mode  = 1;
param.use_gpu        = 0;
param.eigreg         = 1E-4;
%% h-hggm
for k_penalty = 1:length(penalty)
    param.penalty  = penalty(k_penalty);
    [Sjj_est,~,~] = higgs(Svv,Lvj,param);
end

end




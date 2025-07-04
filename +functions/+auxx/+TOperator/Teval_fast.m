function [T,Ginv] = Teval_fast(parameters)
import functions.*
import functions.auxx.*
import functions.auxx.Regularization.*
%% 
Nw = parameters.Dimensions.Nw;
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
if Nv==Nr
    C= nulldiag(parameters.Compact_Model.C);
    K = parameters.Compact_Model.K;
    D = parameters.Compact_Model.D;
else
    C = nulldiag(parameters.Model.C);
    K = parameters.Model.K;
    D = parameters.Model.D;
end
D = 0.5*(D+D');
Kinv = pinv(K);
freq = parameters.Data.freq;
%%

T = zeros(Ne,Nv,Nw);
Ginv = zeros(Nv,Ne,Nw);

if parameters.Parallel.T == 1
    parfor w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       w
       Z = C .* exp(-2 * pi * 1i * freq(w) * D);
       T(:,:,w) = K + K * Z;
       Ginv(:,:,w) = Kinv - Z * Kinv;
    end
else
   for w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       Z = C .* exp(-2 * pi * 1i * freq(w) * D);
       T(:,:,w) = K + K * Z;
       Ginv(:,:,w) = Kinv - Z * Kinv;
    end
end 
end
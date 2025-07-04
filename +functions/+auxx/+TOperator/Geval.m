function [G] = Geval(parameters)
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

Kinv = pinv(K);
freq = parameters.Data.freq;
%%

I = eye(size(C));
G = zeros(Nv,Ne,Nw);
PF = floor(Nw/2)+1;
if parameters.Parallel.T == 1
    parfor w = 1:PF
        %fprintf('-->Frequency %d/%d\n',w,Nw)
        Ts = (I - C .* exp(-2 * pi*1i * freq(w) * D));
        G(:,:,w) = Ts*Kinv;
    end
    parfor w = PF+1:Nw
        %fprintf('-->Frequency %d/%d\n',w,Nw)
        Ts = (I - C .* exp(-2 * pi*1i * freq(w) * D));
        G(:,:,w) = Ts*Kinv;
    end
else
   for w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       Ts = (I - C .* exp(-2 * pi*1i * freq(w) * D));
       G(:,:,w) = Ts*Kinv;
    end
end 

end
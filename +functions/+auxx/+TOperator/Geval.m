function G = Geval(parameters);

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


freq = parameters.Data.freq;
%%

I = eye(size(C));
G = zeros(Nv,Nv,Nw);

if parameters.Parallel.T == 1
    parfor w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       Ts = (I - C .* exp(-2 * pi*i * freq(w) * D));
       G(:,:,w) = inv(Ts);
    end
else
   for w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       Ts = (I - C .* exp(-2 * pi*i * freq(w) * D));
       G(:,:,w) = inv(Ts);
    end
end 

end
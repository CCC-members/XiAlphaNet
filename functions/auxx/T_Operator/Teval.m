function  [T] = Teval(parameters)
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
T = zeros(Ne,Nv,Nw);
tic
if parameters.Parallel.T == 1
    parfor w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       Ts = (I - C .* exp(-2 * pi*i * freq(w) * D));
       T(:,:,w) = (K/Ts);
    end
else
   for w = 1:Nw
       %fprintf('-->Frequency %d/%d\n',w,Nw)
       Ts = (I - C .* exp(-2 * pi*i * freq(w) * D));
       T(:,:,w) = (K/Ts);
    end
end 
toc

end
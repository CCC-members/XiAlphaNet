function  [T] = Teval(parameters)
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

rho_C = max(abs(eig(C))); % spectral radius of C
freq = parameters.Data.freq;
%%

I = eye(size(C));
T = zeros(Ne,Nv,Nw);
PF = floor(Nw/2);
if Nv>Nr
    if parameters.Parallel.T == 1
        parfor w = 1:PF
            %fprintf('-->Frequency %d/%d\n',w,Nw)
            Z = C .* exp(-2 * pi * 1i * freq(w) * D);
            if rho_C < 0.1
                T(:,:,w) = K + K * Z;
            else
                T(:,:,w) = K/(I-Z);
            end
        end
        parfor w = 1+PF:Nw
            %fprintf('-->Frequency %d/%d\n',w,Nw)
            Z = C .* exp(-2 * pi * 1i * freq(w) * D);
            if rho_C < 0.1
                T(:,:,w) = K + K * Z;
            else
                T(:,:,w) = K/(I-Z);
            end
        end
    else
        for w = 1:Nw
            %fprintf('-->Frequency %d/%d\n',w,Nw)
            Z = C .* exp(-2 * pi * 1i * freq(w) * D);
            if rho_C < 0.1
                T(:,:,w) = K + K * Z;
            else
                T(:,:,w) = K/(I-Z);
            end
        end
    end
else
    for w = 1:Nw
        %fprintf('-->Frequency %d/%d\n',w,Nw)
        Z = I-C .* exp(-2 * pi * 1i * freq(w) * D);
        T(:,:,w) = K / Z;
    end
end
end
function [Sjj]=mn_cross(Svv,K,Kclass)
[~,~,Nw] = size(Svv);
[N1,N2,Nk] = size(K);
%K = reshape(K,[Ne,Nw,Nk]);
if Kclass == 0
    Sjj = zeros(N2,N2,Nw); % Compute the Cross Spectra
    if Nk == 1
        U_map = pinv(K);
    else
        for i = 1:Nw
            U_map(:,:,i) = pinv(K(:,:,i));
        end
    end
    for i=1:Nw
        if Nk ==1
            Sjj(:,:,i) = U_map*Svv(:,:,i)*U_map';
        else
            Sjj(:,:,i) = U_map(:,:,i)*Svv(:,:,i)*U_map(:,:,i)';
        end
    end
else
    Sjj = zeros(max(N1,N2),Nw); % Compute only the power spectra
    if Nk == 1
        for i=1:Nw
            Svi = Svv(:,:,i);             % Nv × Nv
            Svi = 0.5*(Svi+Svi')';
            M   = K * Svi;               % Ne × Nv
            Sjj(:,i) = real(diag( M * K'));
        end
    else
        for i = 1:Nw
            Ki  = K(:,:,i);               % Ne × Nv
            Svi = Svv(:,:,i);             % Nv × Nv
            M   = Ki * Svi;               % Ne × Nv
            Sjj(:,i) = real(diag(M*Ki'));%real(sum( M .* Ki, 2 ));
        end
    end

end

end

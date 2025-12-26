function I_F = evaluateFisher(x,Ne,Nv,T,sw,sp,nsf_band,Sw)

import functions.*
import functions.auxx.*
import functions.auxx.ModelVectorization.*
import functions.auxx.XiAlphaGradient.*
import functions.auxx.OptimizedOperations.*

% Unpack parameters
[e,a,sigma2] = x2v(x);

Nsw = length(sw(1,:));
I = eye(Ne);

% Parameter sizes
Ne_xi    = 3*Nv;
Ne_alpha = 4*Nv;
Ne_sig   = 1;

Npar = Ne_xi + Ne_alpha + Ne_sig;
I_F  = zeros(Npar,Npar);

% Storage
dXi_de     = zeros(Nv,3,Nsw);
dAlpha_da = zeros(Nv,4,Nsw);
Mdiag     = zeros(Nv,Nsw);
traceR2   = zeros(1,Nsw);

%% Loop over frequencies
for j = 1:Nsw

    omega = sw(1,j);
    p     = sp(j);
    wj    = sw(2,j)/nsf_band;

    S_omega = Sw(:,:,p);
    T_omega = T(:,:,p);

    % Spectral generators
    xi_omega = e(:,1) ./ (1 + e(:,2).*omega.^2).^e(:,3);
    alpha_omega = 0.5 * ( ...
        a(:,1) ./ (1 + a(:,2).*(omega - a(:,4)).^2).^a(:,3) + ...
        a(:,1) ./ (1 + a(:,2).*((-omega) - a(:,4)).^2).^a(:,3) );

    % Derivatives
    [dXi_de(:,:,j), dAlpha_da(:,:,j)] = dXiAlpha_deda(e,a,omega);

    % Covariance
    R_omega = sigma2*I + computeTDT(T_omega,xi_omega + alpha_omega);
    R_omega = regularize_tensor(R_omega);

    invR = pinv(R_omega);

    % Core contraction matrix
    M = T_omega' * invR * T_omega;
    Mdiag(:,j) = real(diag(M));

    % For noise block
    traceR2(j) = real(trace(invR*invR));

end

%% Assemble Fisher matrix

% --- Xi–Xi block ---
for i = 1:Nv
    for k = 1:3
        idx_i = (i-1)*3 + k;
        for j = 1:Nv
            for l = 1:3
                idx_j = (j-1)*3 + l;
                val = 0;
                for w = 1:Nsw
                    val = val + ...
                        Mdiag(i,w)^2 * ...
                        dXi_de(i,k,w) * dXi_de(j,l,w) * ...
                        sw(2,w)/nsf_band;
                end
                I_F(idx_i,idx_j) = val;
            end
        end
    end
end

% --- Alpha–Alpha block ---
offset_a = Ne_xi;
for i = 1:Nv
    for k = 1:4
        idx_i = offset_a + (i-1)*4 + k;
        for j = 1:Nv
            for l = 1:4
                idx_j = offset_a + (j-1)*4 + l;
                val = 0;
                for w = 1:Nsw
                    val = val + ...
                        Mdiag(i,w)^2 * ...
                        dAlpha_da(i,k,w) * dAlpha_da(j,l,w) * ...
                        sw(2,w)/nsf_band;
                end
                I_F(idx_i,idx_j) = val;
            end
        end
    end
end

% --- Xi–Alpha block ---
for i = 1:Nv
    for k = 1:3
        idx_i = (i-1)*3 + k;
        for j = 1:Nv
            for l = 1:4
                idx_j = offset_a + (j-1)*4 + l;
                val = 0;
                for w = 1:Nsw
                    val = val + ...
                        Mdiag(i,w)^2 * ...
                        dXi_de(i,k,w) * dAlpha_da(j,l,w) * ...
                        sw(2,w)/nsf_band;
                end
                I_F(idx_i,idx_j) = val;
                I_F(idx_j,idx_i) = val;
            end
        end
    end
end

% --- Noise block ---
I_F(end,end) = sum(traceR2 .* sw(2,:)/nsf_band);

end

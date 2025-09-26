function [T, G] = Teval(parameters)
import functions.*
import functions.auxx.*
import functions.auxx.Regularization.*

Nw = parameters.Dimensions.Nw;
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;

% --- Choose compact or full model ---
if Nv == Nr
    %% === Case 1: Compact model – always double ===
    C = nulldiag(parameters.Compact_Model.C);
    K = parameters.Compact_Model.K;
    D = parameters.Compact_Model.D;
    freq = parameters.Data.freq;
    I = eye(size(C));
    T = zeros(Ne, Nv, Nw);
    G = zeros(Nv, Ne, Nw);
    Kinv = pinv(K);

    if parameters.Parallel.T == 1
        parfor w = 1:Nw
            Z = C .* exp(-2*pi*1i*freq(w) * D);
            Ts = (I - Z);
            T(:,:,w) = K / Ts;
            G(:,:,w) = Ts * Kinv;
        end
    else
        for w = 1:Nw
            Z = C .* exp(-2*pi*1i*freq(w) * D);
            Ts = (I - Z);
            T(:,:,w) = K / Ts;
            G(:,:,w) = Ts * Kinv;
        end
    end

else
    %% === Case 2: Large model – try single, fallback to double ===
    C = single(nulldiag(parameters.Model.C));
    K = single(parameters.Model.K);
    D = single(parameters.Model.D);
    freq = single(parameters.Data.freq);
    I = eye(size(C), 'single');
    T = zeros(Ne, Nv, Nw, 'single');
    G = zeros(Nv, Ne, Nw, 'single');
    Kinv = pinv(K);

    % Wrap constants for parfor
    Cconst = parallel.pool.Constant(@() C);
    Dconst = parallel.pool.Constant(@() D);

    % Stability check
    Z0 = C .* exp(-2*pi*1i*freq(1) * D);
    rc_global = rcond(double(I - Z0));

    if rc_global > 1e-6
        % --- Single precision safe ---
        if parameters.Parallel.T == 1
            parfor w = 1:Nw
                C = Cconst.Value; D = Dconst.Value;
                Z = C .* exp(-2*pi*1i*freq(w) * D);
                Ts = (I - Z);
                T(:,:,w) = K / Ts;
                G(:,:,w) = Ts * Kinv;
            end
        else
            for w = 1:Nw
                C = Cconst.Value; D = Dconst.Value;
                Z = C .* exp(-2*pi*1i*freq(w) * D);
                Ts = (I - Z);
                T(:,:,w) = K / Ts;
                G(:,:,w) = Ts * Kinv;
            end
        end

    else
        % --- Ill-conditioned: fallback to double ---
        warning('System ill-conditioned (rcond=%.2e). Recomputing in double.', rc_global);
        delete(Cconst); delete(Dconst);

        C = nulldiag(parameters.Model.C);
        K = parameters.Model.K;
        D = parameters.Model.D;
        freq = parameters.Data.freq;
        I = eye(size(C));
        T = zeros(Ne, Nv, Nw);
        G = zeros(Nv, Ne, Nw);
        Kinv = pinv(K);

        PF = floor(Nw/2);

        if parameters.Parallel.T == 1
            parfor w = 1:PF
                Z = C .* exp(-2*pi*1i*freq(w) * D);
                Ts = (I - Z);
                T(:,:,w) = K / Ts;
                G(:,:,w) = Ts * Kinv;
            end
            parfor w = PF+1:Nw
                Z = C .* exp(-2*pi*1i*freq(w) * D);
                Ts = (I - Z);
                T(:,:,w) = K / Ts;
                G(:,:,w) = Ts * Kinv;
            end
        else
            for w = 1:Nw
                Z = C .* exp(-2*pi*1i*freq(w) * D);
                Ts = (I - Z);
                T(:,:,w) = K / Ts;
                G(:,:,w) = Ts * Kinv;
            end
        end
    end
end

%% === FORCE DOUBLE OUTPUT ===
T = double(T);
G = double(G);

end

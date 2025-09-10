function x = generateRandomSample_fit(Nv, Svv, K, freq, var,parameters)
    import functions.*
    import functions.auxx.*
    import functions.auxx.GenerateSourceSample.*
    import functions.auxx.ModelVectorization.*
    import functions.auxx.Simulations.*
    % Compute cross-spectrum
    [Sjj] = mn_cross(Svv, K,1);

    % Preallocate
    ma = zeros(Nv, 4);  % Alpha parameters
    me = zeros(Nv, 3);  % Xi parameters
    Nw = length(freq);
    Sjj_full = zeros(Nv, Nw);

    % Use identity rotation if Nv == 360


    % Diagonalize Sjj across frequencies
    [~,~,Nk] = size(Sjj);
    if Nk==1
        Sjj_full = Sjj;
    else
        parfor j = 1:Nw
            Sjj_full(:, j) = real(diag(Sjj(:, :, j)));
        end
    end


    % Fit Xi-Alpha model per source
    parfor j = 1:Nv
        %fprintf("Fitting ROI %d/%d...\n", j, Nv);
        try
            parms_j = fit_xi_alpha_multi(Sjj_full(j, :)', freq(:), var, 0);  % doplot = 0
            me(j, :) = parms_j(1:3);   % Xi parameters
            ma(j, :) = parms_j(4:7);   % Alpha parameters
        catch ME
            warning("Fit failed at ROI %d: %s", j, ME.message);
            me(j, :) = NaN;
            ma(j, :) = NaN;
        end
    end

    % Final conversion to model vector x
    s2 = 1;
    x = v2x(me, ma, s2);
end





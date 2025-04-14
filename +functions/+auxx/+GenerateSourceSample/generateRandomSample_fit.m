function x = generateRandomSample_fit(Nv,Nr, Svv, K, R, freq, var)
    import functions.*
    import functions.auxx.*
    import functions.auxx.GenerateSourceSample.*
    import functions.auxx.ModelVectorization.*
    import functions.auxx.Simulations.*
    % Compute cross-spectrum
    [Sjj] = mn_cross(Svv, K);

    % Preallocate
    ma = zeros(Nr, 4);  % Alpha parameters
    me = zeros(Nr, 3);  % Xi parameters
    Nw = length(freq);
    Sjj_full = zeros(Nr, Nw);

    % Use identity rotation if Nv == 360
    

    % Diagonalize Sjj across frequencies
    parfor j = 1:Nw
            Sjj_full(:, j) = real(diag(Sjj(:, :, j)));
    end

    % Fit Xi-Alpha model per source
    parfor j = 1:Nr
        %fprintf("Fitting ROI %d/%d...\n", j, Nr);
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
    if Nv>Nr
        me0(:,1) =R'*me(:,1); 
        me0(:,1) = mean(me(:,1)).* me0(:,1)/max(me0(:,1));
        me0(:,2) =R'*me(:,2); 
        me0(:,2) = mean(me(:,2)).* me0(:,2)/max(me0(:,2));
        me0(:,3) =R'*me(:,3); 
        me0(:,3) = mean(me(:,3)).* me0(:,3)/max(me0(:,3));
        ma0(:,1) =R'*ma(:,1); 
        ma0(:,1) = mean(ma(:,1)).*ma0(:,1)/max(ma0(:,1));
        ma0(:,2) =R'*ma(:,2); 
        ma0(:,2) = mean(ma(:,2)).* ma0(:,2)/max(ma0(:,2));
        ma0(:,3) =R'*ma(:,3); 
        ma0(:,3) = mean(ma(:,3)).* ma0(:,3)/max(ma0(:,3));
        ma0(:,4) =R'*ma(:,4); 
        ma0(:,4) = mean(ma(:,4)).* ones(size(R,2),1);
    else 
        me0 =me;
        ma0 =ma;
    end

    x = v2x(me0, ma0, s2);
end

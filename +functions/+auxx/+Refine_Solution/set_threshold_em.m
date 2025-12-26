function threshold = set_threshold_em(V)

    % -------------------------------
    % Basic sanitation
    % -------------------------------
    V = V(~isnan(V));
    V = V(:);
    V = V(abs(V) > 1e-10);

    % Hard fallback if too few points
    if numel(V) < 10
        threshold = prctile(V, 90);
        return
    end

    % Hard fallback if nearly constant
    if std(V) < 1e-6
        threshold = max(V);
        return
    end

    % -------------------------------
    % Robust scaling (improves EM)
    % -------------------------------
    Vm = median(V);
    Vs = mad(V,1);
    if Vs > 0
        Vz = (V - Vm) / Vs;
    else
        threshold = prctile(V, 90);
        return
    end

    % -------------------------------
    % Attempt 2-component GMM
    % -------------------------------
    try
        gm = fitgmdist( ...
            Vz, 2, ...
            'Replicates', 10, ...
            'CovarianceType', 'diagonal', ...
            'SharedCovariance', false, ...
            'RegularizationValue', 1e-2, ...
            'Options', statset('MaxIter', 1000, 'TolFun', 1e-4) ...
        );
    catch
        threshold = prctile(V, 90);
        return
    end

    % -------------------------------
    % Extract and sort components
    % -------------------------------
    mu = gm.mu;
    sigma = squeeze(gm.Sigma);
    w = gm.ComponentProportion;

    % Reject degenerate mixture
    if any(w < 0.05) || any(sigma < 1e-4)
        threshold = prctile(V, 90);
        return
    end

    [mu, idx] = sort(mu);
    sigma = sigma(idx);
    w = w(idx);

    % -------------------------------
    % Gaussian intersection
    % -------------------------------
    pdf1 = @(x) w(1) * exp(-(x-mu(1)).^2 ./ (2*sigma(1))) ./ sqrt(2*pi*sigma(1));
    pdf2 = @(x) w(2) * exp(-(x-mu(2)).^2 ./ (2*sigma(2))) ./ sqrt(2*pi*sigma(2));

    f = @(x) pdf1(x) - pdf2(x);

    a = mu(1) - 5*sqrt(sigma(1));
    b = mu(2) + 5*sqrt(sigma(2));

    if f(a)*f(b) < 0
        z_th = fzero(f, [a b]);
        threshold = z_th * Vs + Vm;   % back to original scale
    else
        threshold = prctile(V, 90);
    end

end

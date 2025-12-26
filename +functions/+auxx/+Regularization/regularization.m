function Xreg = regularization(X)
% REGULARIZATION
% PSD regularization with loop (d < 1000) and closed-form (d >= 1000)

% ------------------------------------------------------------
% Enforce symmetry
% ------------------------------------------------------------
X = (X + X.') / 2;
[d, ~] = size(X);

% ------------------------------------------------------------
% Common quantities
% ------------------------------------------------------------
L    = trace(X) / d;
Xreg = X;

% ============================================================
% SMALL MATRICES: loop-based regularization
% ============================================================
if d < 1000

    emin  = min(eig(X));
    p     = 0;
    delta = 1e-4;

    while emin < 0 && p < 1
        Xreg = (1 - p) * X + p * L * eye(d);
        emin = min(eig(Xreg));
        p    = p + delta;
    end

% ============================================================
% LARGE MATRICES: closed-form + iterative eigen solver
% ============================================================
else

    opts.tol   = 1e-6;
    opts.maxit = 500;

    try
        lambda_min = eigs(X, 1, 'SA', opts);   % smallest algebraic
    catch
        lambda_min = eigs(X, 1, 'SA');
    end

    if lambda_min < 0
        denom = L - lambda_min;
        if denom > eps
            p = -lambda_min / denom;
            Xreg = (1 - p) * X + p * L * speye(d);
        end
    end

end

end

function lambda = lambda_regspace2(freq, T, Sw, index_stoch, Nsfreq, x0)
% Compute regularization parameter space lambda values

% --- Imports 
import functions.auxx.StochasticEval.*
import functions.auxx.GenerateSourceSample.*
import functions.FunctionGrandientProx.*
import functions.auxx.RegSpace.*
import functions.auxx.ModelVectorization.*

[Ne, Nv, ~] = size(T);

% Sample frequencies
[nsf_band, sw, sp] = sample_frequencies(freq, index_stoch, Nsfreq);
epsilon = 10^(-6);
% Initialize parameters
[e, a, s2] = x2v(x0);
e(:,1)     = 1;       % fix e
e(:,2:end) = epsilon;
a(:,1)     = 1;       % fix a
a(:,2:end) = epsilon;
s2         = epsilon;       % variance off

x_init = v2x(e, a, s2);

% Evaluate gradient at this scale
dF = evaluatedF(x_init, Ne, Nv, T, sw, sp, nsf_band, Sw);
[dFe, dFa, dFs2] = x2v(dF);

% Compute lambdas for each block
lambda_blocks = [ ...
    compute_lambda(dFs2);  % Ridge
    compute_lambda(dFe);   % Group Lasso (e)
    compute_lambda(dFa)    % Group Lasso (a)
];

% Apply safety margin to avoid trivial solution
lambda = 0.9 * (lambda_blocks);

end

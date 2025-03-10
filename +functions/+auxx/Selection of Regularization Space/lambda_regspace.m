function  [lambda] = lambda_regspace(freq,T,Sw,L,index_stoch,Nsfreq)
[Ne,Nv,~] = size(T);
[nsf_band, sw, sp] = sample_frequencies(freq, index_stoch, Nsfreq);
x0  = generateRandomSample(Ne, Nv, Sw, T, freq, 0.01);
epsilon = 10^(-1);
dim = size(x0);
x0 = ones(dim)*epsilon;
dF = evaluatedF(x0, Ne,Nv, T, sw, sp, nsf_band, Sw);
[dFe,dFa,dFs2] = x2v(dF);
lambda(1) = compute_lambda(dFs2,L);
lambda(2) = compute_lambda(dFe,L);
lambda(3) = compute_lambda(dFa,L);
lambda = lambda*9/10; % select the 10% of the regression space to avoid the null solution
end

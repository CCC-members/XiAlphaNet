function  [lambda] = lambda_regspace(parameters)
L = parameters.Lipschitz;
x0  = generateRandomSample(parameters,0.001);
dim = size(x0);
epsilon = 10^(-1);
dF = evaluatedF(ones(dim)*epsilon,parameters);
[dFe,dFa,dFs2] = x2v(dF);
lambda(1) = compute_lambda(dFs2,L);
lambda(2) = compute_lambda(dFe,L);
lambda(3) = compute_lambda(dFa,L);
lambda = lambda*9/10; % select the 10% of the regression space to avoid the null solution
end

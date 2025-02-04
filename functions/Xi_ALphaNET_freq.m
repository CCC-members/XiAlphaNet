function [x] =  Xi_ALphaNET_freq(parameters)

disp('-->> Initializing Model Parameters...')
parameters.Parallel.Lambda = 0;
parameters.Parallel.Kmin = 0;
parameters.Parallel.StochFISTA = 0;
parameters.Parallel.Lipt = 1;

disp('-->> Estimating Number of Batchs to StochFISTA...')
k_min = 30;
parameters.Parallel.Kmin = 0;
parameters.Stochastic.stoch = 1;
parameters.Stochastic.Nsfreq = k_min;
parameters.Stochastic.Niter = 1;
parameters = sample_frequencies(parameters);
parameters.Threshold = activation_threshold(parameters);

disp('-->> Estimating Lipschitz Constant...')
parameters.Parallel.Lipt = 1;
parameters.Lipschitz = 0.01;%estimateLipschitzConstant(parameters, 1, 20);
lambda_space  = [100,1000,1000];

disp('-->> Initializing Bayesian Optimization On Regularization...')
parameters.BayesIter = 100;
parameters.Parallel.Lambda = 1;
[lambda_opt] = bayesianOptSearch(lambda_space, parameters);
parameters.Stochastic.stoch = 0;
parameters.Stochastic.Niter = 1;
parameters.Lipschitz=parameters.Lipschitz;
parameters.Threshold=parameters.Threshold;
parameters = sample_frequencies(parameters);

disp('-->> Initializing Stochastic FISTA global optimizer...')
parameters.Parallel.StochFISTA = 0;
[x_opt,~] = stoch_fista_global(lambda_opt, parameters);
x.Solution = x_opt.Solution;
x.Lambda_reg = lambda_opt;
x.K_min= k_min;
end

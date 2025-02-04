function [x,temp_parameters] =  Xi_ALphaNET(parameters)

disp('-->> Initializing Model Parameters...')
%Nworkers = 30;
%parpool(Nworkers);

parameters.Dimensions.Nv = 360;
parameters.BayesIter0 = 100;
parameters.BayesIter = 20;
parameters.Parallel.Conn_Delay = 1;
parameters.Parallel.Lambda = 0;
parameters.Parallel.T = 0;
parameters.Parallel.Kmin = 0;
parameters.Parallel.StochFISTA = 0;
parameters.Parallel.Lipt = 1;
%std_lambdad = 0.57;% The ratio between the standar deviation of the conduction delays over the mean across the lifespan
lambda_space_cd = [[0.43,1.57];[0.001,1.2]]; 

disp('-->> Estimating Connectivity & Delays Strengths...')
[lambda_opt_dc] = bayes_search_conn_delay(lambda_space_cd, parameters);

lambda1 = lambda_opt_dc(1); % Estimated delay strenght
lambda2 = lambda_opt_dc(2); % Estimated connectivity delay
% 
% % Use the lambda values in your calculations
% %parameters = parameters; % Create a copy to avoid modifying the original
% parameters.Model.D = lambda1 * parameters.Model.D;
% parameters.Model.C = lambda2 * parameters.Model.C;
% 
% parameters = Teval(parameters);
% 
% disp('-->> Estimating Number of Batchs to StochFISTA...')
% k_min = 30;
% parameters.Parallel.Kmin = 0;
% parameters.Stochastic.stoch = 1;
% parameters.Stochastic.Nsfreq = k_min;
% parameters.Stochastic.Niter = 1;
% parameters = sample_frequencies(parameters);
% parameters.Threshold = activation_threshold(parameters);
% 
% disp('-->> Estimating Lipschitz Constant...')
% parameters.Parallel.Lipt = 1;
% parameters.Lipschitz = 0.01;%estimateLipschitzConstant(parameters, 1, 100);
% lambda_space  = [100,1000,1000];
% 
% disp('-->> Initializing Bayesian Optimization On Regularization...')
% parameters.BayesIter = 100;
% parameters.Parallel.Lambda = 1;
% [lambda_opt] = bayesianOptSearch(lambda_space, parameters);
% parameters.Stochastic.stoch = 0;
% parameters.Stochastic.Niter = 1;
% parameters.Lipschitz=parameters.Lipschitz;
% parameters.Threshold=parameters.Threshold;
% parameters = sample_frequencies(parameters);
% 
% disp('-->> Estimating Transfer Function...')
% 
% parameters.Dimensions.Nv = 8003;
% parameters.Parallel.T = 1;
% parameters = Teval(parameters);
% temp_parameters = parameters;
% 
% disp('-->> Initializing Stochastic FISTA global optimizer...')
% parameters.Parallel.StochFISTA = 0;
% [x_opt,~] = stoch_fista_global(lambda_opt, parameters);
% x.Solution = x_opt.Solution;
x.Lambda_DC = lambda_opt_dc;
x.Lambda_reg = lambda_opt;
x.K_min= k_min;
end
% 
% % % 
% figure(4)
% x0 = x.Solution;
% s = reconstructSDM(x0,temp_parameters);
% [l0] = log_spectrum(s,temp_parameters);
% plot(temp_parameters.Data.freq,l0');
% figure(5) 
% [l] = log_spectrum(temp_parameters.Data.Cross,temp_parameters);
% plot(temp_parameters.Data.freq,l');
% norm(l0-l,'fro')/norm(l,'fro')
function [x] =  Xi_ALphaNET_delays(parameters)
disp('-->> Initializing Model Parameters...')
parameters.Dimensions.Nv = 360;
parameters.BayesIter0 = 10;
parameters.BayesIter = 10;
parameters.Parallel.Conn_Delay = 0;
parameters.Parallel.Lambda = 0;
parameters.Parallel.T = 0;
parameters.Parallel.Kmin = 0;
parameters.Parallel.StochFISTA = 0;
parameters.Parallel.Lipt = 0;
lambda_space_cd = [[0.4,1.8];[0.001,2]]; 
disp('-->> Estimating Connectivity & Delays Strengths...')
[lambda_opt_dc] = bayes_search_conn_delay(lambda_space_cd, parameters);
x.Lambda_DC = lambda_opt_dc;
end

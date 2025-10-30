%% Main 
% clear all
%load('Data/Model_Parameters/parameters.mat');
parameters.Dimensions.Nv = 8003;
parameters.BayesIter =10;
parameters.Parallel.Conn_Delay = 0;
parameters.Parallel.Lambda = 0;
parameters.Parallel.T = 0;
parameters.Parallel.Kmin = 1;
parameters.Parallel.StochFISTA = 0;
parameters.Parallel.Lipt = 1;

k_min = findMinimumK(parameters,10,20);
parameters.BayesIter=30;
parameters.Stochastic.stoch = 1;
parameters.Stochastic.Nsfreq =k_min;
parameters.Stochastic.Niter = 1;
parameters = sample_frequencies(parameters);
parameters.Threshold = activation_threshold(parameters);
parameters.Lipschitz = 400;%estimateLipschitzConstant(parameters,2,20); 
lambda_space = [100,1000,1000];%lambda_regspace(paters);%

%%
%stoch_fista_global(lambda, Ne,Nv,T,freq,index_stoch,index_parll,Nsfreq,Cross,Nrand,Lipschitz)
lambda= [100,1000,1000];
Ne = 18;
Nv = 8003;
T = parameters.Model.T;
index_stoch = 1;
index_parll = 1;
Nsfreq = 18;
Cross = parameters.Data.Cross;
Nrand = 10;
Lipschitz = 400;

disp('-->> Initializing Stochastic FISTA  global optimazer...')
[lambda_opt] = [100,1000,1000];%bayesianOptSearch(lambda_space, parameters);
parameters.Stochastic.stoch = 1; 

[x_opt,hist] = stoch_fista_global(lambda, Ne,Nv,T,freq,index_stoch,index_parll,Nsfreq,Cross,Nrand,Lipschitz);

x_opt = np_ref_solution(x_opt);
x = x_opt.Solution;

[e,a,s2] = x2v(x);

% Create a single figure for all subplots
%figure;

% % Alpha Activation Plot
% subplot(2, 3, 1); % 3 rows, 1 column, 1st subplot
% J = a(:,1);
% esi_plot(gca,J); % Pass the current axes handle to esi_plot
% title('Alpha Activation');
% 
% % Xi Activation Plot
% subplot(2, 3, 2); % 3 rows, 1 column, 2nd subplot
% J = e(:,1);
% esi_plot(gca,J); % Pass the current axes handle to esi_plot
% title('Xi Activation');
% 
% % Alpha Peak Frequency
% subplot(2, 3, 3); % 3 rows, 1 column, 3rd subplot
% J = a(:,4).*a(:,1)/max(a(:,1));
% esi_plot(gca,J); % Pass the current axes handle to esi_plot
% title('Alpha Peak Freq');
%%
%parameters.Stochastic.stoch = 0;
%parameters = sample_frequencies(parameters);
%[cxi,calpha] = eval_source_con(x,parameters);

%%

 S= reconstructSDM(x, parameters);
 figure(1)
 log_spec = log_spectrum(S,parameters);
  plot((log_spec)')
  figure(2)
  log_spec = log_spectrum(parameters.Data.Cross,parameters);
   plot((log_spec)')
% 
% figure(2)
% [J]=source_cross_dummy(parameters.Data.Cross,parameters);
% log_spec = log_spectrum(J,parameters);
% plot(mean(log_spec)')


%% Main 
% clear all
%load('Data/Model_Parameters/parameters.mat');
parameters.Dimensions.Nv = 8003;
parameters.BayesIter =10;
parameters.Parallel.Conn_Delay = 0;
parameters.Parallel.Lambda = 1;
parameters.Parallel.T = 0;
parameters.Parallel.Kmin = 1;
parameters.Parallel.StochFISTA = 0;
parameters.Parallel.Lipt = 1;

%k_min = findMinimumK(parameters,10,20);
parameters.BayesIter=30;
parameters.Stochastic.stoch = 1;
parameters.Stochastic.Nsfreq =25;
parameters.Stochastic.Niter = 1;
parameters = sample_frequencies(parameters);
parameters.Threshold = activation_threshold(parameters);
parameters.Lipschitz = 0.01;%estimateLipschitzConstant(parameters,2,20); 
lambda_space = [100,1000,1000];%lambda_regspace(paters);%
% tic;
disp('-->> Initializing Stochastic FISTA  global optimazer...')
[lambda_opt] = bayesianOptSearch(lambda_space, parameters);
parameters.Stochastic.stoch = 0;  
tic
[x_opt,hist] = stoch_fista_global(lambda_opt,parameters);
toc
x_opt = np_ref_solution(x_opt);
x = x_opt.Solution;
%toc;
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
 plot(log_spec')
 
figure(2)
[J]=source_cross_dummy(parameters.Data.Cross,parameters);
log_spec = log_spectrum(J,parameters);
plot(parameters.Data.freq,mean(log_spec)')


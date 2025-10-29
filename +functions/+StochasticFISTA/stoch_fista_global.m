function [x_opt, History] = stoch_fista_global(lambda, Ne,Nv,T,freq,index_stoch,index_parll,Nsfreq,Cross,Nrand,Lipschitz,xx0)

import functions.*
import functions.FunctionGrandientProx.*
import functions.auxx.*
import functions.auxx.StochasticEval.*
import functions.auxx.GenerateSourceSample.*
import functions.StochasticFISTA.*

% Extract number of stochastic iterations
% Nw = parameters.Dimensions.Nw;
% Ne = parameters.Dimensions.Ne;
% Nv = parameters.Dimensions.Nv;
% T = parameters.Model.T;
% sw = parameters.Stochastic.Sampled.sw;    % stoc sampled value of freq
% sp = parameters.Stochastic.Sampled.sp;           % stoch sampled position of freq
% nsf_band = parameters.Stochastic.Nsfreq; % number of freq stoch sampled in a band
% Sw = parameters.Data.Cross;


%% Get T



%Nrand = parameters.Stochastic.Niter;
[nsf_band,sw,sp] = sample_frequencies(freq,index_stoch,Nsfreq);
max_backtracking_iter = 60;
L0 = Lipschitz;
max_iter = 100;
tol = 10^(-2);
eta = 2;
var = 0.001; % Variance
f = @(x) evaluateF(x,Ne,T,sw,sp,nsf_band,Cross);
g1 = @(x) evalg(x,1);
g2 = @(x) evalg(x,2);
g3 = @(x) evalg(x,3);
g = @(x) [g1(x), g2(x), g3(x)]';
grad_f = @(x) evaluatedF(x,Ne,Nv,T,sw,sp,nsf_band,Cross);
prox2 = @(x, c) prox(x, c);
History = cell(1, Nrand);

% To store results from parallel execution
x_opts = cell(1, Nrand);
F_vals = inf(1, Nrand);

% Iterate over each random simulation in parallel
if index_parll == 1
    parfor j = 1:Nrand
        tol = 10^(-2);
        x0 = generateRandomSample(xx0, var);%Ne,Nv,Cross,T,freq, var);
        [x, hist, FSmooth] = fista_with_backtracking(f, grad_f, g, prox2, x0, lambda, L0, eta, max_iter, tol, max_backtracking_iter);
        History{j} = hist;

        % Store results of each iteration
        x_opts{j}.Solution = x;
        x_opts{j}.Initguess = x0;
        x_opts{j}.Feval = hist(end); % Corrected to use the last value in 'hist'
        x_opts{j}.FSmooth = FSmooth;
        F_vals(j) = hist(end);

        % Display current simulation status
        %fprintf('Simulation %d/%d, OF Value: %f, Smooth OF Value: %e\n', j, Nrand, hist(end), FSmooth);
    end
else
    for j = 1:Nrand
        x0 = generateRandomSample(xx0, var);
        [x, hist, FSmooth] = fista_with_backtracking(f, grad_f, g, prox2, x0, lambda, L0, eta, max_iter, tol, max_backtracking_iter);
        History{j} = hist;

        % Store results of each iteration
        x_opts{j}.Solution = x;
        x_opts{j}.Initguess = x0;
        x_opts{j}.Feval = hist(end); % Corrected to use the last value in 'hist'
        x_opts{j}.FSmooth = FSmooth;
        F_vals(j) = hist(end);

        % Display current simulation status
       % fprintf('Simulation %d/%d, OF Value: %f, Smooth OF Value: %e\n', j, Nrand, hist(end), FSmooth);
    end
end
% After the parallel loop, find the best solution
[~, best_idx] = min(F_vals);
x_opt = x_opts{best_idx};
%best_idx

end



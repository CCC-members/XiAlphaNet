function [x_opt, History] = stoch_fista_global(lambda, parameters)
    % Extract number of stochastic iterations
    Nrand = parameters.Stochastic.Niter;
    parameters = sample_frequencies(parameters);
    max_backtracking_iter = 60;
    L0 = parameters.Lipschitz;
    max_iter = 60;
    tol = 1e-2;
    eta = 2;
    var = 2; % Variance
    f = @(x) evaluateF(x, parameters);
    g1 = @(x) evalg1(x, parameters);
    g2 = @(x) evalg2(x, parameters);
    g3 = @(x) evalg3(x, parameters);
    g = @(x) [g1(x), g2(x), g3(x)]';
    grad_f = @(x) evaluatedF(x, parameters);
    prox2 = @(x, c) prox(x, c,parameters);
    History = cell(1, Nrand);
    F = inf;
   
    % To store results from parallel execution
    x_opts = cell(1, Nrand); 
    F_vals = inf(1, Nrand); 
   tic
    % Iterate over each random simulation in parallel
    if parameters.Parallel.StochFISTA == 1
        parfor j = 1:Nrand
            x0 = generateRandomSample(parameters, var);
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
            x0 = generateRandomSample(parameters, var);
            [x, hist, FSmooth] = fista_with_backtracking(f, grad_f, g, prox2, x0, lambda, L0, eta, max_iter, tol, max_backtracking_iter);
            History{j} = hist;
    
            % Store results of each iteration
            x_opts{j}.Solution = x;
            x_opts{j}.Initguess = x0;
            x_opts{j}.Feval = hist(end); % Corrected to use the last value in 'hist'
            x_opts{j}.FSmooth = FSmooth;
            F_vals(j) = hist(end);
    
            % Display current simulation status
            fprintf('Simulation %d/%d, OF Value: %f, Smooth OF Value: %e\n', j, Nrand, hist(end), FSmooth);
        end
    end
    % After the parallel loop, find the best solution
    [~, best_idx] = min(F_vals);
    x_opt = x_opts{best_idx};
toc
end


    
% Define the objective function components and their gradient/proximal operators
x0 = zeros(2, 1);
A = rand(2);
b = rand(2,1)
f = @(x) 0.5 * norm(A * x - b)^2;
grad_f = @(x) A' * (A * x - b);
g = @(x) norm(x, 1); % L1 regularization (LASSO)
prox = @(x, tau) sign(x) .* max(0, abs(x) - tau); % Proximal operator for L1 norm

% Parameters
lambda = 0.1;
L0 = 1;
eta = 2;
max_iter = 1000;
tol = 1e-6;

% Run FISTA with backtracking
[x, history] = fista_with_backtracking(f, grad_f, g, prox, x0, lambda, L0, eta, max_iter, tol,300);

% Plot the objective function value over iterations
figure;
plot(history);
xlabel('Iteration');
ylabel('Objective function value');
title('FISTA with Backtracking');

function [x, history, smoothf] = fista_with_backtracking(f, grad_f, g, prox, x0, lambda, L0, eta, max_iter, tol, max_backtracking_iter)
    % f: Smooth part of the objective function
    % grad_f: Gradient of the smooth part of the objective function
    % g: Nonsmooth part of the objective function (regularization term)
    % prox: Proximal operator of the regularization term
    % x0: Initial point
    % lambda: Regularization parameter
    % L0: Initial value of L
    % eta: Backtracking parameter
    % max_iter: Maximum number of iterations
    % tol: Tolerance for stopping criterion based on the function value change
    % max_backtracking_iter: Maximum number of backtracking iterations
    
    % Initialize variables
    x = x0;
    y = x0;
    t = 1;
    L = L0;
    history = zeros(max_iter, 1); % To store the function value at each iteration
    f_old = f(x) + lambda * g(x); % Initial function value

    for k = 1:max_iter
        % Backtracking line search
        found_L = false;
        backtracking_iter = 0;
        fy = f(y);
        gy = g(y);
        grad_fy = grad_f(y);
        while ~found_L && backtracking_iter < max_backtracking_iter
            L_bar = L * eta^backtracking_iter;
            x_new = prox(y - (1 / L_bar) * grad_fy, lambda / L_bar);
            fx_new = f(x_new);
            gx_new = g(x_new);
            Q_L = fy + grad_fy' * (x_new - y) + (L_bar / 2) * norm(x_new - y)^2 + lambda * gx_new;
            if fx_new + lambda * gx_new <= Q_L
                found_L = true;
                L = L_bar;
            else
                backtracking_iter = backtracking_iter + 1;
            end
        end
        
        % Update x, t, and y based on the restart mechanism
        f_new = fx_new + lambda * gx_new;
        if f_new > f_old  % Check if the function value increased
            % Restart mechanism
            t = 1;
            y = x;  % Restart y to the last best x
        else
            % Regular FISTA update
            t_old = t;
            t = (1 + sqrt(1 + 4 * t^2)) / 2;
            y = x_new + ((t_old - 1) / t) * (x_new - x);
            x = x_new; % Update x to new x only if f_new <= f_old
        end

        % Store the function value at the current iteration
        history(k) = f_new;
        
        % Check stopping criterion based on relative change in function value
        if abs(f_new - f_old) / max(1, abs(f_old)) < tol
            break;
        end
        f_old = f_new; % Update old function value for next iteration
    end
    smoothf = f(x);
    history = history(1:k); % Trim the history to the number of iterations performed
end

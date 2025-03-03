function testComputeDiagonalElements
    % Number of simulations
    k = 100;  % Change this to adjust the number of simulations

    % Test parameters
    n = 19; % Number of rows in T
    m = 8000; % Number of columns in T and rows in D
    k_D = 8000; % Number of columns in D

    % Initialize error storage
    errors_non_diagonal = zeros(k, 1);
    errors_diagonal = zeros(k, 1);

    % Perform k simulations
    for sim = 1:k
        % Generate random matrices T and D
        T = rand(n, m);
        D = rand(m, k_D);

        % Test with non-diagonal D
        computed_diag_elements = computeDiagonalElements(T, D, 0);
        expected_diag_elements = diag(T * D * T');
        errors_non_diagonal(sim) = norm(computed_diag_elements - expected_diag_elements, 'fro');

        % Generate a diagonal matrix D
        D_diag = diag(rand(m, 1));

        % Test with diagonal D
        
        disp('-->Optimize time')
        computed_diag_elements_diag = computeDiagonalElements(T, D_diag, 1);
        
        disp('-->Brute time')
        expected_diag_elements_diag = diag(T * D_diag * T');
        
        errors_diagonal(sim) = norm(computed_diag_elements_diag - expected_diag_elements_diag, 'fro');
    end

    % Display average errors
    fprintf('Average error for Non-Diagonal D over %d simulations: %f\n', k, mean(errors_non_diagonal));
    fprintf('Average error for Diagonal D over %d simulations: %f\n', k, mean(errors_diagonal));
end



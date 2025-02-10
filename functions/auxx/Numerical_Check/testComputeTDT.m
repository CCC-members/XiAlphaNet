function testComputeTDT
    % Number of simulations
    k = 100;

    % Test parameters
    n = 19; % Number of rows in T
    m = 8000; % Number of columns in T, assuming square T for simplicity

    % Initialize error storage
    errors = zeros(k, 1);

    % Perform k simulations
    for sim = 1:k
        % Generate random matrix T
        T = rand(n, m);

        % Generate a random diagonal matrix D represented as a vector
        D = diag(rand(1, m));

        % Compute TDT' using the custom function
        tic;
        disp('--> Optimal time')        
        TDT = computeTDT(T, D);
        toc;

        % Compute TDT' using direct matrix multiplication for verification
        tic;
        disp('-->Brute Time')
        TDT_direct = T * D* T';
        toc;

        % Compute the error as the Frobenius norm of the difference
        errors(sim) = norm(TDT - TDT_direct, 'fro');
    end

    % Display the mean error
    mean_error = mean(errors);
    fprintf('Mean error over %d simulations: %f\n', k, mean_error);
end
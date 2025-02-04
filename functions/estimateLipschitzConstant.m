function L_est = estimateLipschitzConstant(parameters, var, numSamples)
    % Estimate the Lipschitz constant of the gradient of a function F
    %
    % Inputs:
    %   parameters - parameters used in evaluateF and generateRandomSample
    %   var - variance parameter for sample generation
    %   numSamples - number of random samples to generate for estimation
    %
    % Outputs:
    %   L_est - estimated Lipschitz constant

    maxRatio = 0;  % Initialize the maximum ratio
    tic
    if parameters.Parallel.Lipt == 0
        for i = 1:numSamples
            %fprintf('-->Lipschitz step %d\n',i)
            % Generate two random samples
            x = generateRandomSample(parameters, var);
            y = generateRandomSample(parameters, var);
            
            % Evaluate the gradients at x and y
            [dFx, Fx, smoothFx] = evaluatedF(x, parameters);
            [dFy, Fy, smoothFy] = evaluatedF(y, parameters);
            
            % Calculate the norms of the gradient difference and the point difference
             gradDiffNorm = norm(dFx - dFy,'fro');
             pointDiffNorm = norm(x - y,'fro');
        %    gradDiffNorm = max(abs(dFx - dFy));
         %   pointDiffNorm = max(abs(x - y));
    
            % Calculate the ratio of the norms
            if pointDiffNorm > 0  % Avoid division by zero
                currentRatio = gradDiffNorm / pointDiffNorm;
                % Update the maximum ratio found
                maxRatio = max(maxRatio, currentRatio);
            end
        end
    else
        parfor i = 1:numSamples
        %fprintf('-->Lipschitz step %d\n',i)
        % Generate two random samples
        x = generateRandomSample(parameters, var);
        y = generateRandomSample(parameters, var);
        
        % Evaluate the gradients at x and y
        [dFx, Fx, smoothFx] = evaluatedF(x, parameters);
        [dFy, Fy, smoothFy] = evaluatedF(y, parameters);
        
        % Calculate the norms of the gradient difference and the point difference
         gradDiffNorm = norm(dFx - dFy,'fro');
         pointDiffNorm = norm(x - y,'fro');
    %    gradDiffNorm = max(abs(dFx - dFy));
     %   pointDiffNorm = max(abs(x - y));

        % Calculate the ratio of the norms
        if pointDiffNorm > 0  % Avoid division by zero
            currentRatio = gradDiffNorm / pointDiffNorm;
            % Update the maximum ratio found
            maxRatio = max(maxRatio, currentRatio);
        end
        end
    end
    toc
    % The estimated Lipschitz constant is the maximum ratio found
    L_est = maxRatio;
    fprintf('-->Lipschitz constant found %d\n',L_est)
end

function L_est = estimateLipschitzConstant(freq,T,Sw,index_stoch,Nsfreq,index_parall, var, numSamples,x0)
import functions.auxx.StochasticEval.*
import functions.auxx.GenerateSourceSample.*
import functions.FunctionGrandientProx.*
    % Estimate the Lipschitz constant of the gradient of a function F
    %
    % Inputs:
    %   parameters - parameters used in evaluateF and generateRandomSample
    %   var - variance parameter for sample generation
    %   numSamples - number of random samples to generate for estimation
    %
    % Outputs:
    %   L_est - estimated Lipschitz constant
    [Ne,Nv,~] = size(T);
    [nsf_band, sw, sp] = sample_frequencies(freq, index_stoch, Nsfreq);
    maxRatio = 0;  % Initialize the maximum ratio
    
    if index_parall == 0
        for i = 1:numSamples
            %fprintf('-->Lipschitz step %d\n',i)
            % Generate two random samples
            x = generateRandomSample(x0, var); 
            y = generateRandomSample(x0, var); 
            
            % Evaluate the gradients at x and y
            [dFx, Fx, smoothFx] = evaluatedF(x, Ne,Nv, T, sw, sp, nsf_band, Sw);
            [dFy, Fy, smoothFy] = evaluatedF(y, Ne,Nv, T, sw, sp, nsf_band, Sw);
            
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
            x = generateRandomSample(x0, var); 
            y = generateRandomSample(x0, var); 
            
            % Evaluate the gradients at x and y
            [dFx, Fx, smoothFx] = evaluatedF(x, Ne,Nv, T, sw, sp, nsf_band, Sw);
            [dFy, Fy, smoothFy] = evaluatedF(y, Ne,Nv, T, sw, sp, nsf_band, Sw);
            
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
    
    % The estimated Lipschitz constant is the maximum ratio found
    L_est = maxRatio;
    fprintf('-->Lipschitz constant found %d\n',L_est)
end

function [eParams,aParams] = fitXiAlpha_eachChannel(parameters)
    % fitXiAlpha_eachChannel
    % Fits the Xi-Alpha model to each channel in the dataset using lsqcurvefit.

    % Extract data
    S = parameters.Data.Cross;
    [N_channels, ~, N_frequencies] = size(S);
    freq = parameters.Data.freq(:);  % Ensure freq is a column vector

    % Initialize N_all by extracting diagonals
    N_all = zeros(N_channels, N_frequencies);
    for j = 1:N_frequencies
        S_j = S(:, :, j);
        if size(S_j, 1) ~= N_channels || size(S_j, 2) ~= N_channels
            error(['S(:,:,', num2str(j), ') is not square or does not match N_channels.']);
        end
        N_all(:, j) = real(diag(S_j));
    end

    % Initialize parameter storage
    eParams = zeros(N_channels, 3);
    aParams = zeros(N_channels, 4);

    % Define bounds
    lb = [0, 0, 0,  0, 0, 0, 0];         % all >= 0
    ub = [Inf, Inf, Inf, Inf, Inf, Inf, max(freq)];

    % Initial guess
    paramInit = [4, 0.01, 1.5, 8, 0.04, 2.0, 45];

    % Optimization options
    options = optimset('Display','off');

    % Loop over each channel
    for iCh = 1:N_channels
        try
            % Extract measured data
            N_measured = N_all(iCh, :)';  % Column vector (N_frequencies x 1)

            % Define the model handle
            fun = @(p, x) xiAlphaModel_oneChannel(p, x);

            % Debugging: Check sizes
            model_output_init = xiAlphaModel_oneChannel(paramInit, freq);
            if ~isequal(size(model_output_init), size(N_measured))
                error(['Size mismatch in channel ', num2str(iCh), ': ', ...
                       'Model output size ', mat2str(size(model_output_init)), ...
                       ' does not match N_measured size ', mat2str(size(N_measured))]);
            end

            % Perform the fit
            [paramFit, resnorm] = lsqcurvefit(fun, paramInit, freq, N_measured, lb, ub,options);

            % Store fitted parameters
            eParams(iCh, :) = paramFit(1:3);
            aParams(iCh, :) = paramFit(4:7);

            % Display summary
%             fprintf('\nChannel %d fitted parameters:\n', iCh);
%             fprintf('  e1 = %.4f, e2 = %.4f, e3 = %.4f\n', eParams(iCh,1), eParams(iCh,2), eParams(iCh,3));
%             fprintf('  a1 = %.4f, a2 = %.4f, a3 = %.4f, a4 = %.4f\n', ...
%                      aParams(iCh,1), aParams(iCh,2), aParams(iCh,3), aParams(iCh,4));
%             fprintf('  Residual norm = %.6f\n', resnorm);

%             % (Optional) Plot measured vs. fitted
%             figure; hold on; grid on;
%             plot(freq, N_measured, 'bo', 'DisplayName','Data');
%             plot(freq, xiAlphaModel_oneChannel(paramFit, freq), 'r-', ...
%                  'LineWidth', 1.5, 'DisplayName','Fitted Model');
%             xlabel('Frequency'); ylabel('N(\omega)');
%             legend('Location','best');
%             title(sprintf('Channel %d Fit', iCh));

        catch ME
            fprintf('Error fitting channel %d: %s\n', iCh, ME.message);
            eParams(iCh, :) = NaN;
            aParams(iCh, :) = NaN;
            
            % (Optional) Handle plotting for failed fits
            % figure; hold on; grid on;
            % plot(freq, N_measured, 'bo', 'DisplayName','Data');
            % xlabel('Frequency'); ylabel('N(\omega)');
            % legend('Location','best');
            % title(sprintf('Channel %d Fit Failed', iCh));
        end


    end

%     % Display all fitted parameters
%     disp('All channels fitted. Parameter matrices:');
%     disp('eParams = '); disp(eParams);
%     disp('aParams = '); disp(aParams);
end

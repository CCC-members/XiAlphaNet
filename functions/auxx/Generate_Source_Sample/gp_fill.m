function [delay] = gp_fill(delay_mat,length_mat);

% reaging problem size

[N,~] = size(delay_mat);

% First extract the non nan data into a vector;

[delay_vec,length_vec]  = extract_non_nan(delay_mat,length_mat);

%Finding the zeros 
zeros_vec = delay_vec.*length_vec;

non_zero_index = (zeros_vec~=0);

delay_vec = delay_vec(non_zero_index);
length_vec = length_vec(non_zero_index);

% Map the vector into log space 

log_delay_vec = log(delay_vec);
log_length_vec = log(length_vec);

% Gaussian Process regression 

% Set up optimization options
%opts = struct('Optimizer', 'bayesopt', 'ShowPlots', true, 'Verbose', 1, ...
  %            'AcquisitionFunctionName', 'expected-improvement-plus');
gprMdl = fitrgp(log_length_vec, log_delay_vec, 'Basis', 'constant', 'FitMethod', 'exact',...
                'PredictMethod', 'exact', 'KernelFunction', 'squaredexponential');

%plot
xPred = linspace(min(log_length_vec), max(log_length_vec), 100)';  % Points at which to make predictions
[yPred, yStd] = predict(gprMdl, xPred);  % Predictions and standard deviations

% Plotting
figure(4);
errorbar(xPred, yPred, 1.96*yStd, 'k');  % 95% confidence interval
hold on;
index = randperm(200);
plot(log_length_vec(index), log_delay_vec(index), 'ro');  % Original data
title('Gaussian Process Regression with Optimized Hyperparameters');
xlabel('Tract Length log(L_{ij})');
ylabel('Tract Delay log(D_{ij})');
legend('95% Prediction Interval', 'Data', 'Location', 'Best');

for i=1:N
    for j=1:N
        if isnan(delay_mat(i,j))
            delay_mat(i,j) = exp(predict(gprMdl, log(length_mat(i,j))));
        end
    end
end
delay = delay_mat;
end

function ComputeTensorField(min_cond_delay)
%% Tensor_Field
% Read initial parameters
clc;clear all;
load('Data/Model_Parameters/parameters.mat');
parameters.Model.D=0.0110*parameters.Model.D/ mean(parameters.Model.D(:));
output_dir = fullfile('Data','Tensor_Field2');
parameters.Parallel.T = 1;
% Define Grid of Delays and Connectivity
lambda_D_min = 0.4;
lambda_D_max = 1.6;
lambda_C_min = 0.01;
lambda_C_max = 2;
% Number of point of the grid
num_D = 11;
num_C = 11;
% Define Grid
lambda_D_space = linspace(lambda_D_min,lambda_D_max,num_D);
lambda_C_space = linspace(lambda_C_min,lambda_C_max,num_C);
% Generate Tensor Field
for i = 1:num_D
    i
    for j=1:num_C
        tic
        % Extract current point in the field 
        lambda_D = lambda_D_space(i);
        lambda_C = lambda_C_space(j);
        current_parameters = parameters;
        % Scale Delays and Connectivity;
        current_parameters.Model.D = lambda_D*parameters.Model.D;
        current_parameters.Model.C = lambda_C*parameters.Model.C;
        % Evaluate the Tensor Field 
        current_parameters = Teval(current_parameters);
        T = current_parameters.Model.T;
        % Define the filename using the grid indices
        filename = sprintf('T_%d_%d',i,j);
        filepath = fullfile(output_dir,filename);
        save(filepath,'T','-v7.3');
        clear  T
        clear current_parameters
        toc
    end
end

end


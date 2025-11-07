function ComputeTensorField(age,parameters,output_dir)
%% Tensor_Field
import app.*
import app.functions.*
import functions.*
import functions.StochasticFISTA.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.CrossSpectrum.*
import functions.auxx.ModelVectorization.*
import functions.auxx.Regularization.*
import functions.auxx.TOperator.*
import functions.auxx.GenerateSourceSample.*
import functions.auxx.RegSpace.*
import functions.auxx.Simulations.*
import tools.*

if age < 15
parameters.Model.D=0.0110*parameters.Model.D/ mean(parameters.Model.D(:));
end
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
    for j=1:num_C 
        % Extract current point in the field 
        lambda_D = lambda_D_space(i);
        lambda_C = lambda_C_space(j);
        current_parameters = parameters;
        % Scale Delays and Connectivity;
        current_parameters.Model.D = lambda_D*parameters.Model.D;
        current_parameters.Model.C = lambda_C*parameters.Model.C;
        % Evaluate the Tensor Field 
        T = Teval(current_parameters);
        % Define the filename using the grid indices
        filename = sprintf('T_%d_%d',i,j);
        filepath = fullfile(output_dir,filename);
        save(filepath,'T','-v7.3');
        clear  T
        clear current_parameters
    end
end

end


function loss = delay_objective_function(delay,age_star,h,batch,parameters)
%% Compute the T operator for current delay
parameters.Model.D = v2m(delay);
parameters =  Teval(parameters);
%% 
file_path_data ='Data/Scalp_Density_Matrix/Control/';
matFiles = dir(fullfile(file_path_data, '*.mat'));
rand_list = randperm(length(matFiles),batch);
loss=0;
for j =1:length(rand_list)
    filePath = fullfile(file_path_data, matFiles(rand_list(j)).name);
    %filePath_x = fullfile(folderPath_x, matFiles_x(j).name);

    % Load the .mat file
    data_struct = load(filePath);
    
    %% Extract and store the required data and parameters
    parameters.Data.Cross = data_struct.data_struct.CrossM(:,:,1:49);
    parameters.Data.Age = data_struct.data_struct.age;
    if ischar(data_struct.data_struct.age) || isstring(data_struct.data_struct.age)
        % Convert the string to a number
        age = str2double(data_struct.data_struct.age);
    else
        % If it's already a number, just use it directly
        age = data_struct.data_struct.age;
    end
    parameters.Data.Age = age;
    parameters.Stochastic.stoch = 1;
    parameters.Stochastic.Nsfreq =1;
    parameters.Stochastic.Niter = 1;
    parameters = sample_frequencies(parameters);
    parameters.Threshold = activation_threshold(parameters);
    parameters.Lipschitz = 0.01;%estimateLipschitzConstant(parameters, 1, 60);
    lambda_space = [100,1000,1000];%lambda_regspace(parameters);
    disp('-->> Initializing Stochastic FISTA  global optimazer...')
    %[lambda_opt] = bayesianOptSearch(lambda_space, parameters);
    lambda_opt = [4.8,248.26,689.26];
    parameters.Stochastic.stoch = 1;
    parameters.Stochastic.Nsfreq =1;
    parameters.Stochastic.Niter = 1;
    parameters = sample_frequencies(parameters);
    [x_opt,~] = stoch_fista_global(lambda_opt,parameters);
    x = x_opt.Solution;
    loss  = loss + Kernel_age(age_star,age,h)*evaluateF(x,parameters);
end
loss = loss/batch;
end


function [T,pos] = read_transfer_function_tensor_field(delay,conn,age,TF_path)
%load("Data/Tensor_Field/");
% Define Grid of Delays and Connectivity
lambda_D_min = 0.4;
lambda_D_max = 1.6;
lambda_C_min = 0.01;
lambda_C_max = 2;
% Number of point of the grid
num_D = 10;
num_C = 10;
% Define Grid
lambda_D_space = linspace(lambda_D_min,lambda_D_max,num_D);
lambda_C_space = linspace(lambda_C_min,lambda_C_max,num_C);

%%
[~,pos_delay] = min(abs(delay-lambda_D_space));
[~,pos_conn] = min(abs(conn-lambda_C_space));

pos= [pos_delay,pos_conn];
%%
if age>15
    filename = fullfile(TF_path,'Tensor_Field_9.5ms',strcat("T_",num2str(pos_delay),'_',num2str(pos_conn)));
    T = load(filename);
    T = T.T;
else
    filename = fullfile(TF_path,'Tensor_Field_11ms',strcat("T_",num2str(pos_delay),'_',num2str(pos_conn)));
    T = load(filename);
    T = T.T;
end
end

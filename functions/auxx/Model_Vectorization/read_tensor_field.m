function [T,pos] = read_tensor_field(delay,conn,age)
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
    filename = sprintf("Data/Tensor_Field/T_%d_%d.mat",pos_delay,pos_conn);
    T = load(filename);
    T = T.T;
else
    filename = sprintf("Data/Tensor_Field2/T_%d_%d.mat",pos_delay,pos_conn);
    T = load(filename);
    T = T.T;
end
end

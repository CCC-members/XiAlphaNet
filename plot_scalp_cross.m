TF_path = fullfile(properties.general_params.tmp.path,'TensorField');
figure(4)
x0 = x.Solution;
temp_parameters = parameters;
temp_parameters.Age = data.age;
lambda1 = x.Lambda_DC(1);
lambda2 = x.Lambda_DC(2);
temp_parameters.Model.T = read_tensor_field(lambda1,lambda2,temp_parameters.Age,TF_path);
s = reconstructSDM(x0,temp_parameters,data);
[l0] = log_spectrum(s,data.freq);
plot(data.freq,l0');
figure(5) 
[l] = log_spectrum(data.Cross,data.freq);
plot(data.freq,l');
norm(l0-l,'fro')/norm(l,'fro')
clear temp_parameters
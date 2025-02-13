function [x,T] =  Xi_ALphaNET(properties,data,parameters)

disp('-->> Initializing Model Parameters...')

Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
Ne = parameters.Dimensions.Ne;
Nw = properties.general_params.data.nFreqs;
BayesIter_Delay = properties.model_params.BayesIter_Delay;
BayesIter_Reg1 = properties.model_params.BayesIter_Reg1;
BayesIter_Reg2 = properties.model_params.BayesIter_Reg2;
Nrand1 = properties.model_params.Nrand1;
Nrand2 = properties.model_params.Nrand2;
lambda_space_cd = properties.model_params.delay.lambda_space_cd;
conn_delay = properties.general_params.parallel.conn_delay;
stoch1 = properties.model_params.stoch1;
stoch2 = properties.model_params.stoch2;
tf_default = properties.model_params.tensor_field.default;


parameters.Dimensions.Nv = Nr;
[NewCross,scale] = global_scale_factor_correction(data.Cross);
data.Cross = NewCross;
clear NewCross;

 
if data.age<15
   parameters.Model.D=0.0110*parameters.Model.D/ mean(parameters.Model.D(:));
end
age = data.age;


disp('-->> Estimating Connectivity & Delays Weights...')
Cross = data.Cross;
freq = data.freq;
K = parameters.Compact_Model.K;
D = parameters.Compact_Model.D;
C = parameters.Compact_Model.C;


[lambda_opt_dc] = bayes_search_conn_delay(lambda_space_cd, Ne,Nr,Nw,freq,Cross,BayesIter_Reg1,K,D,C,conn_delay,BayesIter_Delay);


lambda1 = lambda_opt_dc(1); % Estimated delay strenght
lambda2 = lambda_opt_dc(2); % Estimated connectivity delay
% 
% Use the lambda values to updata connectivity and delays

parameters.Model.D = lambda1 * parameters.Model.D;
parameters.Model.C = lambda2 * parameters.Model.C;

parameters.Compact_Model.D = lambda1 * parameters.Compact_Model.D;
parameters.Compact_Model.C = lambda2 * parameters.Compact_Model.C;

parameters.Parallel.T = 0;
parameters.Data.freq = freq;
T = Teval(parameters);
disp('-->> Estimating Number of Batchs to StochFISTA...')

k_min = findMinimumK(freq,T,Cross, 5, 20,conn_delay);

disp('-->> Initializing Bayesian Optimization On Regularization...')
index_parall_bayes= 1;
Nsfreq = k_min;
Lipschitz = 0.01;%estimateLipschitzConstant(freq,T,Cross,1,Nsfreq,stoch1, 0.1, 20);
lambda_space  = [100,1000,1000];%lambda_regspace(freq,T,Cross,Lipschitz,stoch1,Nsfreq);

[lambda_opt] = bayesianOptSearch(lambda_space,Ne,Nr,T,freq,stoch1,0,index_parall_bayes,Nsfreq,Cross,Nrand1,Lipschitz,BayesIter_Reg2);

disp('-->> Estimating Transfer Function...')
if(tf_default)
    TF_path = fullfile(properties.general_params.tmp.path,'TensorField');
    T = read_tensor_field(lambda1,lambda2,age,TF_path);
else
    parameters.Parallel.T = 1;
    T = Teval(parameters);
end
clear parameters;

disp('-->> Initializing Stochastic FISTA global optimizer...')
[x_opt, ~] = stoch_fista_global(lambda_opt, Ne,Nv,T,freq,stoch2,conn_delay,Nsfreq,Cross,Nrand2,Lipschitz);

[e,a,s2] = x2v(x_opt.Solution);
e(:,1) = e(:,1)/scale;
a(:,1) = a(:,1)/scale;
s2 = s2/scale;
x_opt.Solution = v2x(e,a,s2);
x.Solution = x_opt.Solution;
x.Lambda_DC = lambda_opt_dc;
x.Lambda_reg = lambda_opt;
x.Age = data.age;
x.kmin = k_min;
x.Reg_Space = lambda_space;
end

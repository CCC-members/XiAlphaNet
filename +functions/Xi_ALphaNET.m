function [x,T,G] =  Xi_ALphaNET(properties,data,parameters)
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
disp('-->> Initializing Model Parameters...')

Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
Ne = parameters.Dimensions.Ne;
Nw = properties.model_params.nFreqs;
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
 %scale=1;
if data.age<15
   parameters.Model.D=0.0110*parameters.Model.D/ mean(parameters.Model.D(:));
   parameters.Compact_Model.D=0.0110*parameters.Compact_Model.D/ mean(parameters.Compact_Model.D(:));
end
age = data.age;
% 
Cross = data.Cross;
freq = data.freq;
K = parameters.Compact_Model.K;
D = parameters.Compact_Model.D;
C = parameters.Compact_Model.C;
R = parameters.Compact_Model.R;
% 
disp('-->> Fixing Initial parameters...')
T = Teval(parameters);
parameters.Model.T = T;
Cross1 = mn_cross(Cross,K);
for j=1:Nw
    G(:,:,j) = inv(eye(size(C))- C .* exp(-2 * pi* 1i * freq(j) * D));
end
x0 = generateRandomSample_fit(Nr,Nr, Cross1, G, R, freq, 1);

disp('-->> Cross Validating Regularization Parmeters...')
k_min = 30;
index_parall_bayes= 1*conn_delay;
Nsfreq = k_min;
Lipschitz = 10^(14);
% L = 10.^(linspace(0,30,10)); % Space
% Value = [];
% for j = 1:length(L)
%     lambda_space = lambda_regspace(freq,T,Cross,L(j),stoch1,Nsfreq,x0);
%     [lambda_opt] = bayesianOptSearch(lambda_space,Ne,Nr,T,freq,stoch1,1,index_parall_bayes,Nsfreq,Cross,Nrand1,L(j),30,x0);
%     [x_opt, ~] = stoch_fista_global(lambda_opt, Ne,Nv,T,freq,stoch1,1,Nsfreq,Cross,10,L(j),x0);
%     Value(j) = x_opt.Feval;
% end
% [~,index] = min(Value);
% Lipschitz = L(index);


disp('-->> Estimating Connectivity & Delays Weights...')


[lambda_opt_dc] = bayes_search_conn_delay(lambda_space_cd, Ne,Nr,Nw,freq,Cross,BayesIter_Reg1,K,D,C,1,BayesIter_Delay,x0,Lipschitz);


lambda1 = lambda_opt_dc(1); % Estimated delay strenght
lambda2 = lambda_opt_dc(2); % Estimated connectivity delay
% 
% Use the lambda values to updata connectivity and delays

parameters.Model.D = lambda1 * parameters.Model.D;
parameters.Model.C = lambda2 * parameters.Model.C;

parameters.Compact_Model.D = lambda1 * parameters.Compact_Model.D;
parameters.Compact_Model.C = lambda2 * parameters.Compact_Model.C;
C = parameters.Compact_Model.C ;
D = parameters.Compact_Model.D ;


parameters.Parallel.T = 0;
parameters.Data.freq = freq;
T = Teval(parameters);

% disp('-->> Fixing Initial parameters...')
% for j=1:Nw
%     G(:,:,j) = inv(eye(size(C))- C .* exp(-2 * pi*i * freq(j) * D));
% end
% x0 = generateRandomSample_fit(Nr,Nr, Cross1, G, R, freq, 3);


disp('-->> Initializing Bayesian Optimization On Regularization...')
%estimateLipschitzConstant(freq,T,Cross,1,Nsfreq,stoch1, 0.1, 1000,x0);
lambda_space  = lambda_regspace(freq,T,Cross,Lipschitz,stoch1,Nsfreq,x0)

[lambda_opt] = bayesianOptSearch(lambda_space,Ne,Nr,T,freq,stoch1,0,index_parall_bayes,Nsfreq,Cross,Nrand1,Lipschitz,BayesIter_Reg2,x0);

disp('-->> Estimating Transfer Function...')
parameters.Dimensions.Nv = Nv;
if(tf_default)
    TF_path = fullfile(properties.general_params.tmp.path,'TensorField');
    T = read_tensor_field(lambda1,lambda2,age,TF_path);
else
    parameters.Parallel.T = 0;
    T = Teval(parameters);
end
clear parameters;

disp('-->> Initializing Stochastic FISTA global optimizer...')
tic;
[x_opt, ~] = stoch_fista_global(lambda_opt, Ne,Nv,T,freq,stoch2,conn_delay,Nsfreq,Cross,Nrand2,Lipschitz,x0);
toc;
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
function [x, T, G,x0] = Xi_ALphaNET(properties, data, parameters)
% Xi_ALphaNET: Estimates EEG source connectivity, conduction delays, and spectral parameters
%
% INPUTS:
%   properties - structure containing model parameters and general properties
%   data - structure containing EEG cross-spectrum, frequency, and age data
%   parameters - structure containing model dimensions and compact model definitions
%
% OUTPUTS:
%   x - structure with estimated parameters and regularization values
%   T - estimated transfer function tensor
%   G - frequency-dependent connectivity matrix

%% Import required libraries and functions
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

  tic
    disp('-->> Initializing Model Parameters...');
    Nr                          = parameters.Dimensions.Nr;
    Nv                          = parameters.Dimensions.Nv;
    Ne                          = parameters.Dimensions.Ne;
    Nw                          = properties.model_params.nFreqs;
    BayesIter_Delay             = properties.model_params.BayesIter_Delay;
    BayesIter_Reg1              = properties.model_params.BayesIter_Reg1;
    BayesIter_Reg2              = properties.model_params.BayesIter_Reg2;
    Nrand1                      = properties.model_params.Nrand1;
    Nrand2                      = properties.model_params.Nrand2;
    conn_delay                  = properties.general_params.parallel.conn_delay;
    stoch1                      = properties.model_params.stoch1;
    stoch2                      = properties.model_params.stoch2;
    tf_default                  = properties.model_params.tensor_field.default;
    parameters.Dimensions.Nv    = Nr;
    [NewCross,scale]            = global_scale_factor_correction(data.Cross);
    data.Cross                  = NewCross;
    clear NewCross;
    if data.age<15
        parameters.Model.D      = 0.0110*parameters.Model.D/ mean(parameters.Model.D(:));
        parameters.Compact_Model.D  = 0.0110*parameters.Compact_Model.D/ mean(parameters.Compact_Model.D(:));
    end
    age                         = data.age;

    disp('-->> Fixing Initial Parameters...');
    Cross                       = data.Cross;
    freq                        = data.freq;
    K                           = parameters.Compact_Model.K;
    D                           = parameters.Compact_Model.D;
    C                           = parameters.Compact_Model.C;
    parameters.Data.freq        = freq;
    parameters.Parallel.T       = 0;
    T                           = Teval(parameters);
    G                           = Geval(parameters);
    parameters.Model.T          = T;

    R = parameters.Model.R;
    x0                          = generateRandomSample_fit(Nr, Cross, G,G, freq, 10,conn_delay,R);

    disp('-->> Estimating Connectivity & Delays Weights...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda_space_cd             = properties.model_params.delay.lambda_space_cd;
    lambda_space_cd(2,:)        = [10^(-10),1/max(abs(eig(parameters.Model.C)))];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lambda_opt_dc]             = bayes_search_conn_delay(lambda_space_cd, Ne, Nr, Nw, freq, Cross, K, D, C, 1, BayesIter_Delay, x0, 1);
    lambda1                     = lambda_opt_dc(1); % Estimated delay strenght
    lambda2                     = lambda_opt_dc(2); % Estimated connectivity delay

    % Use the lambda values to updata connectivity and delays
    parameters.Model.D          = lambda1 * parameters.Model.D;
    parameters.Model.C          = lambda2 * parameters.Model.C;
    parameters.Compact_Model.D  = lambda1 * parameters.Compact_Model.D;
    parameters.Compact_Model.C  = lambda2 * parameters.Compact_Model.C;

    disp('-->> Estimating Lipschitz Constant...');
    k_min                       = floor(length(freq)*4/5);
    index_parall_bayes          = conn_delay;
    Nsfreq                      = k_min;
    Lipschitz                   = estimateLipschitzConstant(freq, T, Cross, 1, 25, stoch1, 0.001, 100, x0);

    disp('-->> Cross Validating Regularization Space...');
    [lambda_space, ~, ~]        = find_best_lambda(freq, T, Cross, 1, 1, 20, x0, Ne, Nr, Nr, 10, index_parall_bayes, Nrand1, Nrand2, Lipschitz, conn_delay);

    % Initializing Bayesian Optimization On Regularization
    disp('-->> Bayesian Optimization On Regularization Parameters...');
    [lambda_opt] = bayesianOptSearch(lambda_space, Ne, Nr, T, freq, stoch1, 0, index_parall_bayes, Nsfreq, Cross, Nrand1, Lipschitz, BayesIter_Reg2, x0);

    % Estimate transfer function
    disp('-->> Estimating Transfer Function...');
    parameters.Parallel.T       = 1*conn_delay;
    parameters.Dimensions.Nv    = Nv;
    if(tf_default)
        TF_path                 = fullfile(properties.general_params.tmp.path,'TensorField');
        [T,G,~]                 = read_tensor_field(lambda1,lambda2,25,TF_path);
    else
        [T,G]                     = Teval(parameters);
    end

    disp('-->> Fixing Initial Parameters...');

    x0                           = generateRandomSample_fit(Nv, Cross, G,G, freq, 4,conn_delay,R);

    %% Initializing Stochastic FISTA global optimizer
    disp('-->> Running Stochastic FISTA Global Optimization...');
    [x_opt, ~]                  = stoch_fista_global(lambda_opt, Ne, Nv, T, freq, 0, conn_delay, Nsfreq, Cross, Nrand2,Lipschitz, x0);
    [e,a,s2]                    = x2v(x_opt.Solution);
    e(:,1)                      = e(:,1)/scale;
    a(:,1)                      = a(:,1)/scale;
    s2                          = s2/scale;
    data.Cross                  = data.Cross/scale;
    x_opt.Solution              = v2x(e,a,s2);
    x.Solution                  = x_opt.Solution;
    x.Lambda_DC                 = lambda_opt_dc;
    x.Lambda_reg                = lambda_opt;
    x.Age                       = data.age;
    x.kmin                      = k_min;
    x.Reg_Space                 = lambda_space;
    toc

end

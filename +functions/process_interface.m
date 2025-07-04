function  process_interface(properties)

%%
%%  Importing Packages 
%%
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

input_path                              = properties.general_params.input_path;
output_path                             = properties.general_params.output_path;
part_file                               = properties.general_params.dataset.participants_info.file;

%%
%%
%%
XAN_filename                            = 'XIALPHANET.json';
XAN_path                                = fullfile(output_path);
XAN_file                                = fullfile(XAN_path,XAN_filename);
if(isfile(XAN_file))
    XIALPHANET                          = jsondecode(fileread(XAN_file));
    parameters                          = load(fullfile(output_path,XIALPHANET.Structural.parameters));
else
    XIALPHANET.Name                     = properties.general_params.dataset.Name;
    TempUUID                            = java.util.UUID.randomUUID;
    XIALPHANET.UUID                     = char(TempUUID.toString);
    XIALPHANET.Description              = properties.general_params.dataset.Description;
    XIALPHANET.Task                     = properties.general_params.dataset.descriptors.task;
    XIALPHANET.Status                   = "Processing";
    XIALPHANET.Location                 = XAN_path;
    XIALPHANET.general_params           = properties.general_params;    
    XIALPHANET.Participants     = [];

    %%
    %% Preprocessing
    %%
    [XIALPHANET,parameters] = preprocessing(properties, XIALPHANET);
end

%% Loading Participamts Infomation
participants_file = '';
isPartInfo = false;
[~,~,part_ext] = fileparts(fullfile(input_path,part_file));
if(isequal(part_ext,'.tsv') || isequal(part_ext,'.csv'))
    isPartInfo = true;
    opts = detectImportOptions(fullfile(input_path,part_file), FileType="text");
    participants_file = table2struct(readtable(fullfile(input_path,part_file),opts));
end
if(isequal(part_ext,'.json'))
    isPartInfo = true;
    participants_file = jsondecode(fileread(fullfile(input_path,part_file)));
end


%% Loop through each .mat file in the subjects
subjects = dir(input_path);
subjects(ismember({subjects.name},{'..','.','structural'})) = [];
subjects([subjects.isdir]==0) = [];
if(~isempty(properties.general_params.participants))
    subjects = subjects(ismember({subjects.name}, properties.general_params.participants));
end
parameters_tmp = parameters;
lambda_space_cd                 = properties.model_params.delay.lambda_space_cd;
lambda_space_cd(2,:)            = [10^(-10),1/max(abs(eig(parameters.Model.C)))];
for s=1:length(subjects) 
    tic
    if(isfile(XAN_file))
        XIALPHANET              = jsondecode(fileread(XAN_file));
    end
    subject                     = subjects(s);
    SubID                       = subject.name;
    Participant.SubID           = SubID;

    % Update progress, report current estimate
    if(getGlobalGuimode())
        if(properties.dlg.CancelRequested)
            msg = "Are you sure to cancel the processing?";
            title = "Cancel processing";
            answer = uiconfirm(properties.UIFigure,msg,title, ...
                "Options",["Stop","Continue"], ...
                "Icon","+guide/images/question.png",...
                "DefaultOption",1,"CancelOption",2);
            % Handle response
            switch answer
                case 'Stop'
                    break;                    
                case 'Continue'
                    properties.dlg.CancelRequested=false;
                    drawnow;
            end
        end
        properties.dlg.Message = strcat("Processing subject: ",SubID, " (",num2str(s)," of ",num2str(length(subjects)),")");
        drawnow;
    end
      
    %%
    %% Check participant processed
    %%
    if(isempty(XIALPHANET.Participants))
        iPart = 1;
    else
        iPart = find(ismember({XIALPHANET.Participants.SubID},SubID),1);
        if(~isempty(iPart) && isequal(XIALPHANET.Participants(iPart).Status,"Completed"))
            disp('---------------------------------------------------------------------');
            disp(strcat("-->> The analysis of subject: ", SubID, " has been completed."));
            disp(strcat("-->> Jumping to the next subject."));
            disp('---------------------------------------------------------------------');
            continue;
        end
        if(isempty(iPart))
            iPart = length(XIALPHANET.Participants) + 1;
        end
    end
    disp('---------------------------------------------------------------------');
    disp(strcat("-->> Processing subject: ", SubID));
    disp('---------------------------------------------------------------------');

    %%
    %% Check data structure
    %%
    if(isPartInfo)
        if(isfield(participants_file,'participant_id'))
            pat = participants_file(find(ismember({participants_file.participant_id},{SubID}),1));
        end
        if(isfield(participants_file,'SubID'))
            pat = participants_file(find(ismember({participants_file.SubID},{SubID}),1));
        end
    end
    [data,status,Participant] = check_data_structure(properties,Participant,subject,pat);    
    if(~status)
        XIALPHANET.Participants(iPart).SubID = Participant.SubID;
        XIALPHANET.Participants(iPart).Age = Participant.Age;
        XIALPHANET.Participants(iPart).Status = Participant.Status;
        XIALPHANET.Participants(iPart).Errors = Participant.Errors;
        XIALPHANET.Participants(iPart).FileInfo = Participant.FileInfo;
        saveJSON(XIALPHANET,XAN_file);
        disp('---------------------------------------------------------------------');
        continue;
    end

    %%
    %% Initializing Model Parameters
    %
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
    x0                          = generateRandomSample_fit(Nr, Cross, G, freq, 3);

    disp('-->> Estimating Lipschitz Constant...');
    k_min                       = 25;
    index_parall_bayes          = conn_delay;
    Nsfreq                      = k_min;
    Lipschitz                   = estimateLipschitzConstant(freq, T, Cross, 1, 25, stoch1, 0.001, 100, x0);
   
    disp('-->> Cross Validating Initial Regularization Space...');
    [lambda_space, ~, ~]        = find_best_lambda(freq, T, Cross, 1, 1, 20, x0, Ne, Nr, Nr, 10, index_parall_bayes, Nrand1, Nrand2, Lipschitz, conn_delay);

    disp('-->> Estimating Connectivity & Delays Weights...');
    [lambda_opt_dc]             = bayes_search_conn_delay(lambda_space_cd, Ne, Nr, Nw, freq, Cross, BayesIter_Reg1, K, D, C, 1, BayesIter_Delay, x0, Lipschitz, lambda_space);
    lambda1                     = lambda_opt_dc(1); % Estimated delay strenght
    lambda2                     = lambda_opt_dc(2); % Estimated connectivity delay

    % Use the lambda values to updata connectivity and delays
    parameters.Model.D          = lambda1 * parameters.Model.D;
    parameters.Model.C          = lambda2 * parameters.Model.C;
    parameters.Compact_Model.D  = lambda1 * parameters.Compact_Model.D;
    parameters.Compact_Model.C  = lambda2 * parameters.Compact_Model.C;
    parameters.Parallel.T       = 1*conn_delay;
    T                           = Teval(parameters);
    % 
    disp('-->> Cross Validating Final Regularization Space...');
    [lambda_space, ~, ~]        = find_best_lambda(freq, T, Cross, 1, 1, 20, x0, Ne, Nr, Nr, 10, index_parall_bayes, Nrand1, Nrand2, Lipschitz, conn_delay);

    % Initializing Bayesian Optimization On Regularization
    disp('-->> Bayesian Optimization On Regularization Parameters...');
    [lambda_opt] = bayesianOptSearch(lambda_space, Ne, Nr, T, freq, stoch1, 0, index_parall_bayes, Nsfreq, Cross, Nrand1, Lipschitz, BayesIter_Reg2, x0);
    
    % Estimate transfer function
    disp('-->> Estimating Transfer Function...');
    parameters.Dimensions.Nv    = Nv;
    if(tf_default)
        TF_path                 = fullfile(properties.general_params.tmp.path,'TensorField');
        [T]                     = read_tensor_field(lambda1,lambda2,age,TF_path);
        [G]                     = Geval(parameters);
    else
        [T]                     = Teval(parameters);
        [G]                     = Geval(parameters);
    end

    % Initializing Stochastic FISTA global optimizer
    disp('-->> Fixing Initial Parameters...');
    x0                           = generateRandomSample_fit(Nv, Cross, G, freq, 1); 
    disp('-->> Running Stochastic FISTA Global Optimization...');
    [x_opt, ~]                  = stoch_fista_global(lambda_opt, Ne, Nv, T, freq, stoch2, conn_delay, Nsfreq, Cross, Nrand2, Lipschitz, x0);
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

    %%
    %% Saving Participant file
    %%
    [Participant]               = xan_save(properties,'subject',SubID,x,T,G,parameters,data,Participant);
    parameters                  = parameters_tmp;
    toc

    %% Save the computed x to the corresponding group folder in Model_Parameters

    disp('-->> Saving XiAlphaNet Information file.');
    XIALPHANET.Participants(iPart).SubID        = Participant.SubID;
    XIALPHANET.Participants(iPart).Age          = Participant.Age;
    XIALPHANET.Participants(iPart).Status       = "Completed";
    XIALPHANET.Participants(iPart).Errors       = Participant.Errors;
    XIALPHANET.Participants(iPart).FileInfo     = strcat(SubID,".json");
    saveJSON(XIALPHANET,XAN_file);
    disp('---------------------------------------------------------------------');
    disp(strcat("-->> The analysis of subject: ", SubID, " has been completed."));
    disp('---------------------------------------------------------------------');
    if(getGlobalGuimode())
        properties.dlg.Value = s/length(subjects);
        drawnow;
    end
end
%%
%% Group analysis
%%
% XIALPHANET          = groupProcessing(properties,XIALPHANET);

%%
%% Saving XIALPHANET file
%%
XIALPHANET          = jsondecode(fileread(XAN_file));
XIALPHANET.Status   = "Completed";
saveJSON(XIALPHANET,XAN_file);
h                   = matlab.desktop.editor.openDocument(XAN_file);
h.smartIndentContents
h.save
h.close

%%
%%  Updating XiAlphaNet datasets
%%
if(isfile(XAN_file))
    % Loading existed Datasets
    
    Datasets_file = xan_get('datasets_file');
    if(isfile(Datasets_file))
        TempDatasets = jsondecode(fileread(Datasets_file));
        if(~isempty(TempDatasets))
            Datasets = TempDatasets;
        else
            Datasets = struct([]);
        end
    else
        Datasets = struct([]);
    end

    % Including new dataset
    XIALPHANET = jsondecode(fileread(XAN_file));
    if(isempty(Datasets))
        Datasets            = XIALPHANET;
    else
        Datasets(end + 1) = XIALPHANET;
    end
    disp("-->> Dataset saved");
    saveJSON(Datasets,Datasets_file);
    h = matlab.desktop.editor.openDocument(Datasets_file);
    h.smartIndentContents
    h.save
    h.close
end

end

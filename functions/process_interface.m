function  process_interface(properties)

input_path = properties.general_params.input_path;
output_path = properties.general_params.output_path;

%%
%%
%%
XAN_filename = 'XIALPHANET.json';
XAN_path = fullfile(output_path);
XAN_file = fullfile(XAN_path,XAN_filename);
if(isfile(XAN_file))
    XIALPHANET = jsondecode(fileread(XAN_file));
    parameters = load(fullfile(output_path,XIALPHANET.Structural));
else   
    XIALPHANET.Name             = properties.general_params.dataset.Name;    
    TempUUID                    = java.util.UUID.randomUUID;
    XIALPHANET.UUID             = char(TempUUID.toString);
    XIALPHANET.Description      = properties.general_params.dataset.Description;
    XIALPHANET.Task             = properties.general_params.dataset.descriptors.task;
    XIALPHANET.Status           = "Processing";
    XIALPHANET.Location         = XAN_path;
    XIALPHANET.general_params   = properties.general_params;
    XIALPHANET.Structural       = "structural/parameters.mat";
    XIALPHANET.Participants     = [];

    %%
    %% Preprocessing
    %%
    parameters = preprocessing(properties);    
end

%% Loop through each .mat file in the subjects
subjects = dir(input_path);
subjects(ismember({subjects.name},{'..','.','structural'})) = [];
subjects([subjects.isdir]==0) = [];
for s=1:length(subjects)
    if(isfile(XAN_file))
        XIALPHANET              = jsondecode(fileread(XAN_file));
    end
    subject                     = subjects(s);
    SubID                       = subject.name;
    Participant.SubID           = SubID;
    if(isempty(XIALPHANET.Participants))
        iPart = 1;
    else
        iPart = find(ismember({XIALPHANET.Participants.SubID},SubID),1);
        if(~isempty(iPart) && isequal(XIALPHANET.Participants(iPart).Status,"Completed"))
            continue; 
        end
        if(isempty(iPart))
            iPart = length(XIALPHANET.Participants) + 1;
        end
    end
    disp(strcat("-->> Processing subject: ", SubID));
    disp('---------------------------------------------------------------------');
    
    %%
    %% Check data structure
    %%    
    [data,status,Participant] = check_data_structure(properties,Participant,subject);
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
    
    subject_path                    = fullfile(output_path,SubID);
    if(~isfolder(subject_path))
        mkdir(subject_path);
    end

    %%
    %% Analysis level base
    %%


    %%
    %% Analysis level delay
    %%


    %%
    %% Analysis level XXXXXX
    %%
    [x,T] =  Xi_ALphaNET(properties,data,parameters);
   
    %%
    %% Saving Participant file
    %%
    disp('-->> Saving Participant Information file.')
    [e,a,~] = x2v(x.Solution);
    % Saving Modulatory Weights of the Anatomical Connectivity and
    % Conduction Delays
    Mod_Weights = x.Lambda_DC;
    Participant.Mod_Weights = fullfile('Mod_Weights.mat');
    save(fullfile(subject_path,"Mod_Weights.mat"),"Mod_Weights");
    % Saving Parameters of the Xi-Process (Aperiodic Component)
    Xi_estimate.Power = e(:,1);    % Xi - Power
    Xi_estimate.Width = e(:,2);    % Xi - Width
    Xi_estimate.Exponent = e(:,3); % Xi - Exponent
    Participant.Xi_estimate = fullfile('Xi_estimate.mat');
    save(fullfile(subject_path,"Xi_estimate.mat"),'-struct',"Xi_estimate");
    % Saving Parameters of the Alpha-Process 
    Alpha_estimate.Power = a(:,1);    % Alpha - Power
    Alpha_estimate.Width = a(:,2);    % Alpha - Width
    Alpha_estimate.Exponent = a(:,3); % Alpha - Exponent
    Alpha_estimate.PAF = a(:,4);      % Alpha - Peak Frequency
    Participant.Alpha_estimate = fullfile('Alpha_estimate.mat');
    save(fullfile(subject_path,"Alpha_estimate.mat"),'-struct',"Alpha_estimate");
    % Saving Estimated Condunction Delay Matrix between Neurotracts
    Delay_Matrix = x.Lambda_DC(1)*parameters.Compact_Model.D;
    Participant.Delay_Matrix = fullfile('Delay_Matrix.mat');
    save(fullfile(subject_path,"Delay_Matrix.mat"),"Delay_Matrix");
     % Saving Estimated Anatomical Connectivity Matrix between Neurotracts
    Conn_Matrix = x.Lambda_DC(2)*parameters.Compact_Model.C;
    Participant.Conn_Matrix = fullfile('Conn_Matrix.mat');
    save(fullfile(subject_path,"Conn_Matrix.mat"),"Conn_Matrix");
    % Map Solution into Source  Activation and Crosspectrum on the Frequency Domain
    [source_act_cross] = eval_source_conn(x.Solution, data.freq,T,parameters.Model.K,parameters.Model.R,properties); 
    % Saving  Full Activation
    Source_Activations_Full = source_act_cross.Activations.Full;
    Participant.Source_Activations_Full = fullfile('Source_Activations_Full.mat');
    save(fullfile(subject_path,"Source_Activations_Full.mat"),"Source_Activations_Full");
    % Saving Xi-Process Activation
    Source_Activations_Xi = source_act_cross.Activations.Xi;
    Participant.Source_Activations_Xi = fullfile('Source_Activations_Xi.mat');
    save(fullfile(subject_path,"Source_Activations_Xi.mat"),"Source_Activations_Xi");
    % Saving Alpha-Process Activation
    Source_Activations_Alpha = source_act_cross.Activations.Xi;
    Participant.Source_Activations_Alpha = fullfile('Source_Activations_Alpha.mat');
    save(fullfile(subject_path,"Source_Activations_Alpha.mat"),"Source_Activations_Alpha");
    % Saving Full Cross-Spectrum in ROIs-Space
    Source_Cross_Full = source_act_cross.Cross.Full;
    Participant.Source_Cross_Full = fullfile('Source_Cross_Full.mat');
    save(fullfile(subject_path,"Source_Cross_Full.mat"),"Source_Cross_Full");
    % Saving Xi-Process Cross-Spectrum in ROIs-Space
    Source_Cross_Xi = source_act_cross.Cross.Xi;
    Participant.Source_Cross_Xi = fullfile('Source_Cross_Xi.mat');
    save(fullfile(subject_path,"Source_Cross_Xi.mat"),"Source_Cross_Xi");
    % Saving Alpha-Process Cross-Spectrum in ROIs-Space
    Source_Cross_Alpha = source_act_cross.Cross.Alpha;
    Participant.Source_Cross_Alpha = fullfile('Source_Cross_Alpha.mat');
    save(fullfile(subject_path,"Source_Cross_Alpha.mat"),"Source_Cross_Alpha");
    % Saving Status
    Participant.Status = "Completed";
    saveJSON(Participant,fullfile(subject_path ,strcat(SubID,'.json')));

    %% Save the computed x to the corresponding group folder in Model_Parameters
    
    disp('-->> Saving XiAlphaNet Information file.')
    XIALPHANET.Participants(iPart).SubID = Participant.SubID;  
        XIALPHANET.Participants(iPart).Age = Participant.Age;
        XIALPHANET.Participants(iPart).Status = "Completed";
        XIALPHANET.Participants(iPart).Errors = Participant.Errors;
    XIALPHANET.Participants(iPart).FileInfo     = strcat(SubID,".json");
    saveJSON(XIALPHANET,XAN_file);
    disp('---------------------------------------------------------------------');
end


%%
%% Saving XIALPHANET file
%%
XIALPHANET = jsondecode(fileread(XAN_file));
XIALPHANET.Status = "Completed";
saveJSON(XIALPHANET,XAN_file);
h = matlab.desktop.editor.openDocument(XAN_file);
h.smartIndentContents
h.save
h.close


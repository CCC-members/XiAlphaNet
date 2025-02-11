function  process_interface(properties)

input_path = properties.general_params.input_path;
output_path = properties.general_params.output_path;

 
%%
%% Preprocessing
%%
tic;
parameters = preprocessing(properties);
toc;
%%
%%
%%
XAN_filename = 'XIALPHANET.json';
XAN_path = fullfile(output_path);
XAN_file = fullfile(XAN_path,XAN_filename);
if(isfile(XAN_file))
    XIALPHANET = jsondecode(fileread(XAN_file));
else   
    XIALPHANET.Name             = properties.general_params.dataset.Name;    
    TempUUID                    = java.util.UUID.randomUUID;
    XIALPHANET.UUID             = char(TempUUID.toString);
    XIALPHANET.Description      = properties.general_params.dataset.Description;
    XIALPHANET.Task             = properties.general_params.dataset.descriptors.task;
    XIALPHANET.Status           = "Processing";
    XIALPHANET.Location         = XAN_path;
    XIALPHANET.general_params   = properties.general_params;
    XIALPHANET.Participants     = [];
end

%% Loop through each .mat file in the subjects
subjects = dir(input_path);
subjects(ismember({subjects.name},{'..','.'})) = [];
for s=1:10
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
    Participant.Status = "Completed";
    Participant.Estimations = fullfile('x_source_estimations.mat');
    save(fullfile(subject_path,"x_source_estimations.mat"),'-struct',"x");
    Participant.TransferFunction = fullfile('transfer_function.mat');
    save(fullfile(subject_path,'transfer_function.mat'),"T");
    saveJSON(Participant,fullfile(subject_path ,strcat(SubID,'.json')));
    
    % [x] = np_ref_solution(x);
    %[x] = global_scale_factor(x,parameters1);
    
   

    %% Save the computed x to the corresponding group folder in Model_Parameters
    % saveFilePath = fullfile(subject_path,sprintf('x_opt_%d.mat',perml(s)));
    % save(saveFilePath, 'x');
    % clear x
    % clear parameters1

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


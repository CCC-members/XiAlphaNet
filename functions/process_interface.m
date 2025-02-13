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
    [Participant] = xan_save(properties,SubID,'subject',x,T,parameters,Participant);


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


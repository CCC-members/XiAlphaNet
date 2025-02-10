function  process_interface()

properties = get_properties();
input_path = properties.general_params.input_path;
output_path = properties.general_params.output_path;
ref_file = properties.general_params.data.ref_file;
 
%%
%% Preprocessing
%%


%%
%%
%%
% Load previously saved parameters if necessary
load('data/parameters.mat');
%preprocessing_velocity;



%% Loop through each .mat file in the subjects
subjects = dir(input_path);
subjects(ismember({subjects.name},{'..','.'})) = [];
for s=1:length(subjects)

    subject = subjects(s);
    SubID = subject.name;

    disp(strcat("-->> Processing subject: ", SubID));
    data = load(fullfile(subject.folder,subject.name,strrep(ref_file,'SubID',SubID)));
    while isfield(data,'data_struct')
        data = data.data_struct;
    end

    %%
    %% Check data structure
    %%
    [data,status,errors] = check_data_structure(properties,data);
    if(~status)

        continue;
    end

    

    tic
    [x] =  Xi_ALphaNET(properties,data,parameters);
    % [x] = np_ref_solution(x);
    %[x] = global_scale_factor(x,parameters1);
    toc
    x.Age=age;

    %% Save the computed x to the corresponding group folder in Model_Parameters
    saveFilePath = fullfile(output_path, subFolders{k},sprintf('x_opt_%d.mat',perml(s)));
    save(saveFilePath, 'x');
    clear x
    clear parameters1
end


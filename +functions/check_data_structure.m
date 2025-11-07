function [data,status,Participant] = check_data_structure(properties, Participant, subject, varargin)

%%
%%  Importing Packages
%%
import functions.*
import functions.auxx.*
import functions.auxx.DataPreprosessing.*
import functions.auxx.OptimizedOperations.*
import functions.import.*
import plugins.*
import plugins.HarMNqEEG.*

for i=1:length(varargin)
    eval([inputname(i+3) '= varargin{i};']);
end
Nw = properties.model_params.nFreqs;
country = properties.general_params.data.country;
eeg_device = properties.general_params.data.eeg_device;
ref_file = properties.general_params.data.ref_file;
errors = {};
status = true;
type = properties.general_params.data.type;
file_name = fullfile(subject.folder,subject.name,strrep(ref_file,'SubID',Participant.SubID));
try
    switch lower(type)
        case 'cross'
            data = load(file_name);
            while isfield(data,'data_struct')
                data = data.data_struct;
            end
        case 'eeg_signal'
            [data,error_msg] = ImportEEG(properties,file_name,Participant.SubID,pat);
            if(~isempty(error_msg))
                status = false;
                Participant.Age = '';
                Participant.Status = "Rejected";
                Participant.FileInfo = "";
                Participant.Errors{1} = error_msg;
                fprintf(2,strcat('\n-->> Error: The data structure for subject: ',subject.name,' \n'));
                fprintf(2,strcat('-->> Have the folows errors.\n'));
                fprintf(2,strcat("-->> " ,error_msg, ".\n"));
                fprintf(2,strcat('-->> Jump to an other subject.\n'));
                data = [];
                return;
            end
    end

    Participant.Status = "Checked";
    if(~isfield(data,'age'))
        data.age = randi([20,80],1);
    end
    if ischar(data.age) || isstring(data.age)
        % Convert the string to a number
        data.age = str2double(data.age);
    elseif(isnan(data.age))
        data.age = randi([20 80],1);
    else
        data.age = data.age;
    end

    if(~isfield(data,'freqrange'))
        Participant.Status = "Rejected";
        status = false;
    end

    if(~isfield(data,'CrossM'))
        Participant.Status = "Rejected";
        status = false;
    else
        % Cross
        data.Cross = data.CrossM(:,:,1:Nw);
        data.Cross = aveReference(data.Cross);
        data.Cross = regularize_tensor(data.Cross);
        data.freq = data.freqrange(1:Nw);
        data.dnames = data.dnames(1:end-1); % delete reference channel
        data = rmfield(data,{'CrossM','freqrange'});
    end
    Participant.Age = data.age;
    Participant.Errors = errors;
    Participant.FileInfo = "";

catch Ex
    Participant.Age = '';
    Participant.Status = "Rejected";
    Participant.FileInfo = "";
    Participant.Errors{1} = Ex.message;
    status = false;
    fprintf(2,strcat('\n-->> Error: The data structure for subject: ',subject.name,' \n'));
    fprintf(2,strcat('-->> Have the folows errors.\n'));
    for j=1:length(Participant.Errors)
        fprintf(2,strcat("-->> " ,Participant.Errors{j}, ".\n"));
    end
    fprintf(2,strcat('-->> Jump to an other subject.\n'));
    data = [];
    return;
end

if(status)
    disp("-->> Saving subject data");
    SubID = Participant.SubID;
    [Participant] = xan_save(properties,'create_subject',SubID,data,Participant);
else
    fprintf(2,strcat('\n-->> Error: The folder structure for subject: ',subject.name,' \n'));
    fprintf(2,strcat('-->> Have the folows errors.\n'));
    for j=1:length(errors)
        fprintf(2,strcat('-->>' ,errors{j}, '.\n'));
    end
    fprintf(2,strcat('-->> Jump to an other subject.\n'));
end
end


function [data,status,Participant] = check_data_structure(properties, Participant, subject)

Nw = properties.general_params.data.nFreqs;
ref_file = properties.general_params.data.ref_file;
errors = {};
status = true;
try
    data = load(fullfile(subject.folder,subject.name,strrep(ref_file,'SubID',Participant.SubID)));
    while isfield(data,'data_struct')
        data = data.data_struct;
    end
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
    data.freq = data.freqrange(1:Nw); 
    data = rmfield(data,{'CrossM','freqrange'});
end

Participant.Age = data.age;
Participant.Errors = errors;
Participant.FileInfo = "";

if(~status)
    fprintf(2,strcat('\n-->> Error: The folder structure for subject: ',subject.name,' \n'));
    fprintf(2,strcat('-->> Have the folows errors.\n'));
    for j=1:length(errors)
        fprintf(2,strcat('-->>' ,errors{j}, '.\n'));
    end
    fprintf(2,strcat('-->> Jump to an other subject.\n'));
end
end


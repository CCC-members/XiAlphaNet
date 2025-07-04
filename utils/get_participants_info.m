
import tools.*

input_path = "/Users/ronald/Downloads/MultinationalNorms";

subjects = dir(input_path);
subjects(ismember({subjects.name},{'.','..'})) = [];
subjects([subjects.isdir] == 0) = [];

Participants = struct;

for i=1:length(subjects)
    subject = subjects(i);
    disp(strcat("-->> Processing subject: ", subject.name));
    Participants(i).SubID = subject.name;
    try
        load(fullfile(subject.folder,subject.name,strcat(subject.name,'.mat')));
        Participants(i).Age = data_struct.age;
        Participants(i).Sex = data_struct.sex;
    catch
        Participants(i).Age = '';
        Participants(i).Sex = '';
    end
end
saveJSON(Participants,fullfile(input_path,'Participants.json'));


input_path = "";
output_path = "";

files = dir(input_path);
subjects(ismember({files.name},{'..','.'})) = [];

for i=1:length(files)
    file = files(i);
    [~,SubID,~] = fileparts(files);
    disp(strcat("-->> Updating file: ", SubID));
    mkdir(fullfile(output_path,SubID));
    copyfile(fullfile(file.folder,file.name),fullfile(output_path,SubID,file.name));
end
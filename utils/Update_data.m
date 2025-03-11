

input_path = "C:\ronald\last_version_december\New folder\Xi-AlphaNET compress (copy)\Data\Scalp_Density_Matrix\Control";
output_path = "C:\ronald\Data\Norms";

files = dir(input_path);
files(ismember({files.name},{'..','.'})) = [];

for i=1:length(files)
    file = files(i);
    [~,SubID,~] = fileparts(file.name);
    disp(strcat("-->> Updating file: ", SubID));
    mkdir(fullfile(output_path,SubID));
    copyfile(fullfile(file.folder,file.name),fullfile(output_path,SubID,file.name));
end
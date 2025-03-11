input_path = "/mnt/Data/NeuroEPO/Data/Grupo 1 Pre";
output_path = "/mnt/Data/NeuroEPO/Data/Grupo 1 Pre";

files = dir(fullfile(input_path,"*.PLG"));
files(ismember({files.name},{'..','.'})) = [];

for i=1:length(files)
    file = files(i);
    [~,SubID,ext] = fileparts(file.name);
    SubID = strrep(SubID,' ','_');
    disp(strcat("-->> Updating file: ", SubID));
    mkdir(fullfile(output_path,SubID));

    if(isfile(fullfile(file.folder,file.name)))
        movefile(fullfile(file.folder,file.name),fullfile(output_path,SubID,strcat(SubID,'.PLG')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.CDC'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.CDC')),fullfile(output_path,SubID,strcat(SubID,'.CDC')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.CMM'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.CMM')),fullfile(output_path,SubID,strcat(SubID,'.CMM')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.INF'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.INF')),fullfile(output_path,SubID,strcat(SubID,'.INF')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.MRK'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.MRK')),fullfile(output_path,SubID,strcat(SubID,'.MRK')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.pat'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.pat')),fullfile(output_path,SubID,strcat(SubID,'.pat')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.srf'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.srf')),fullfile(output_path,SubID,strcat(SubID,'.srf')),'f');
    end
    if(isfile(fullfile(file.folder,strrep(file.name,ext,'.WIN'))))
        movefile(fullfile(file.folder,strrep(file.name,ext,'.WIN')),fullfile(output_path,SubID,strcat(SubID,'.WIN')),'f');
    end
end
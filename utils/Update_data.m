

input_path = 'D:\XiAlphaNetProject\ParkinsonsXialphaproject\ParkinsonsData\NeuroEPOEEGRawData\Post\Group3Post\';
output_path = 'D:\XiAlphaNetProject\ParkinsonsXialphaproject\ParkinsonsData\DataGatheredNew\Group3Post\';


files = dir(input_path);
files(ismember({files.name},{'..','.'})) = [];

for i=1:length(files)
    file = files(i);
    [~,SubID,~] = fileparts(file.name);
    disp(strcat("-->> Updating file: ", SubID));
    mkdir(fullfile(output_path,SubID));
    copyfile(fullfile(file.folder,file.name),fullfile(output_path,SubID,file.name));
end
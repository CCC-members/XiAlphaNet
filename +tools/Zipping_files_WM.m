

clear all;
close all;
clc;

disp("-->> Starting process.");

% defining vars
addpath(genpath('D:\Tools\cifti-matlab-master'));
addpath(genpath('plugins'));

root_path = "Z:/data3_260T/share_space/jcclab-users/Mitchell/HCP/TimeSeries";
outputpath = "Z:/data3_260T/share_space/jcclab-users/Mitchell/HCP/Zip_WM";

mkdir(outputpath);

subjects = dir(root_path);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

subject_IDs = cell(length(subjects),1);
count = 1;
tic;
load('glasserMaps32K.mat');
parfor i=1:length(subjects)    
    subject = subjects(i);
    subID = subject.name;
    subject_dir = fullfile(root_path,subID);
    files   = dir(fullfile(subject_dir,'tfMRI_WM*'));   
    if(~isempty(files))
        subject_IDs{i} = subID;
        for j=1:length(files)
            file        = files(j);
            old_file    = fullfile(file.folder,file.name); 
            load(old_file, 'dtseries');
            cdata       = detrend(dtseries,'constant');
            gTS         = GlasserTimeSeries(cdata, LGlasserMap, RGlasserMap);
            mkdir(fullfile(outputpath,subID));
            parsave(fullfile(outputpath,subID,file.name),'gTS',gTS);            
        end
        count = count + 1;        
    end
end
save(fullfile(outputpath,'subjects.mat'),'subject_IDs');
zippedfiles = zip(fullfile(outputpath,'fix_WM_gTS.dtseries.zip'),outputpath);
toc;
disp("-->> Process finished.");
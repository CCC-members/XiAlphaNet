

clear all;
close all;
clc;

disp("-->> Starting process.");

% defining vars
addpath(genpath('D:\Tools\cifti-matlab-master'));
addpath(genpath('plugins'));

root_path = "Z:/data3_260T/share_space/jcclab-users/Mitchell/HCP/fix_REST/";
outputpath = "Z:/data3_260T/share_space/jcclab-users/Mitchell/HCP/Zip_fix_REST";

mkdir(outputpath);

subjects = dir(root_path);
subjects(ismember( {subjects.name}, {'.', '..'})) = [];  %remove . and ..

subject_IDs = cell(length(subjects),1);
count = 1;
for i=1:length(subjects)
    tic;
    subject = subjects(i);
    subID = subject.name;
    old_file = fullfile(root_path,subID,'fix_REST_gTS.dtseries.mat');
    if(isfile(old_file))
        subject_IDs{i} = subID;
        mkdir(fullfile(outputpath,subID));
        new_file = fullfile(outputpath,subID,'fix_REST_gTS.dtseries.mat');
        copyfile(old_file,new_file);
        count = count + 1;
    end
end
save(fullfile(outputpath,'subjects.mat'),'subject_IDs');
% zippedfiles = zip('fix_REST_gTS.dtseries.zip',outputpath);

disp("-->> Process finished.");
function xan_create_db()
ciftistormDir = fullfile(getUserDir(),'.CiftiStorm');
if(~isfolder(ciftistormDir))
    mkdir(ciftistormDir);
end
%% Defaults
if(~isfolder(fullfile(ciftistormDir,'defaults')))
    mkdir(fullfile(ciftistormDir,'defaults'));
end
if(~isfolder(fullfile(ciftistormDir,'defaults','anatomy')))
    mkdir(fullfile(ciftistormDir,'defaults','anatomy'));
end
if(~isfolder(fullfile(ciftistormDir,'defaults','eeg')))
    mkdir(fullfile(ciftistormDir,'defaults','eeg'));
    mkdir(fullfile(ciftistormDir,'defaults','eeg','Colin27'));
    mkdir(fullfile(ciftistormDir,'defaults','eeg','ICBM152'));
    mkdir(fullfile(ciftistormDir,'defaults','eeg','NotAligned'));
end
if(~isfolder(fullfile(ciftistormDir,'defaults','meg')))
    mkdir(fullfile(ciftistormDir,'defaults','meg'));
end

%% Datasets
if(~isfolder(fullfile(ciftistormDir,'Datasets')))
    mkdir(fullfile(ciftistormDir,'Datasets'));
end
if(~isfile(fullfile(ciftistormDir,'Datasets','Datasets.json')))
    saveJSON([],fullfile(ciftistormDir,'Datasets','Datasets.json'));
end
end


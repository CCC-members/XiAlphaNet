function  xan_init()
homedir = char(java.lang.System.getProperty('user.home'));
XANdir  = fullfile(homedir,".XiAlphaNet");
if(~isfolder(XANdir))
    mkdir(XANdir);
end
XANDatasetsdir = fullfile(XANdir,'Datasets');
if(~isfolder(XANDatasetsdir))
    mkdir(XANDatasetsdir);
end

XANDatasetsfile = fullfile(XANdir,'Datasets','Datasets.json');
if(~isfile(XANDatasetsfile))
    Datasets = struct([]);
    saveJSON(Datasets,XANDatasetsfile);
end

end


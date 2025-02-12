function xan_set( varargin )

if ((nargin >= 2) && ischar(varargin{1}))
    contextName = varargin{1};
    data = varargin{2};
else
    return
end

% Get required context structure
switch contextName
    case 'cfs_dir'
        
    case 'defaults_dir'
        
    case 'datasets'
        if(~isfolder(fullfile(getUserDir(),'.CiftiStorm')))
            mkdir(fullfile(getUserDir(),'.CiftiStorm'));
        end
        if(~isfolder(fullfile(getUserDir(),'.CiftiStorm','Datasets')))
            mkdir(fullfile(getUserDir(),'.CiftiStorm','Datasets'));
        end
        saveJSON(data,fullfile(getUserDir(),'.CiftiStorm','Datasets','Datasets.json')) 
        h = matlab.desktop.editor.openDocument(fullfile(getUserDir(),'.CiftiStorm','Datasets','Datasets.json'));
        h.smartIndentContents
        h.save
        h.close
    case 'bst_default_eeg'
        
    case 'openmeeg'
        
    case 'duneruro'
end

end
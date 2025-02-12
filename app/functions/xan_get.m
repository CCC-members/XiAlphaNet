function [output1, output2, output3]   = xan_get( varargin )

if ((nargin >= 1) && ischar(varargin{1}))
    contextName = varargin{1};
else
    return
end

% Get required context structure
switch contextName
    case 'cfs_dir'
        output1 = fullfile(getUserDir(),'.CiftiStorm');    
    case 'defaults_dir'
        cfs_db_dir = fullfile(getUserDir(),'.CiftiStorm');
        output1 =  fullfile(cfs_db_dir,'defaults','anatomy');
        output2 =  fullfile(cfs_db_dir,'defaults','eeg');
        output3 =  fullfile(cfs_db_dir,'defaults','meg');
    case 'datasets'
        cfs_db_dir = fullfile(getUserDir(),'.CiftiStorm');        
        datasets_file =  fullfile(cfs_db_dir,'Datasets','Datasets.json');
        if(isfile(datasets_file))
            output1 = jsondecode(fileread(datasets_file));
        else
            output1 = [];
        end
    case 'datasets_file'
        output1 = fullfile(cfs_get( 'cfs_dir' ),'Datasets','Datasets.json');
    case 'bst_default_eeg'
        output1 = 'https://github.com/brainstorm-tools/brainstorm3/raw/master/defaults/eeg';
    case 'openmeeg'
        
    case 'duneruro'
end

end
function [output1, output2, output3]   = xan_get( varargin )

if ((nargin >= 1) && ischar(varargin{1}))
    contextName = varargin{1};
else
    return;
end

% Get required context structure
switch contextName
    case 'xan_dir'
        output1 = fullfile(getUserDir(),'.XiAlphaNet');    
    case 'defaults_dir'
        cfs_db_dir = fullfile(getUserDir(),'.XiAlphaNet');
        output1 =  fullfile(cfs_db_dir,'defaults','anatomy');
        output2 =  fullfile(cfs_db_dir,'defaults','eeg');
        output3 =  fullfile(cfs_db_dir,'defaults','meg');
    case 'datasets'
        cfs_db_dir = fullfile(getUserDir(),'.XiAlphaNet','Datasets');        
        datasets_file =  fullfile(cfs_db_dir,'Datasets','Datasets.json');
        if(isfile(datasets_file))
            output1 = jsondecode(fileread(datasets_file));
        else
            output1 = [];
        end
    case 'datasets_file'
        output1 = fullfile(xan_get( 'xan_dir' ),'Datasets','Datasets.json');
    case 'test_data_url'
        output1 = 'https://github.com/brainstorm-tools/brainstorm3/raw/master/defaults/eeg';
    case 'tensor_field'
        output1 = {'https://github.com/brainstorm-tools/brainstorm3/raw/master/defaults/eeg',
            'https://github.com/brainstorm-tools/brainstorm3/raw/master/defaults/eeg'};
    case 'tmp_path'
        properties = get_properties();
        if(isequal(properties.general_params.tmp.path,'local'))
            output1 = fullfile(pwd,'tmp');
            return;
        end
        if(isfolder(properties.general_params.tmp.path))
            output1 = properties.general_params.tmp.path;
            return;
        end
        output1 = '';
    case 'duneruro'
end

end
function [output1, output2, output3]   = xan_get( varargin )
import app.*
import app.functions.*
import functions.*
import guide.*
import tools.*
if ((nargin >= 1) && ischar(varargin{1}))
    contextName = varargin{1};
    if(nargin >= 2)
        data = varargin{2};
    end
else
    return;
end

% Get required context structure
switch contextName
    case 'xan_dir'
        output1 = fullfile(getUserDir(),'.XiAlphaNet');    
    case 'defaults_dir'
        xan_db_dir = xan_get('xan_db_dir');
        output1 =  fullfile(xan_db_dir,'defaults','anatomy');
        output2 =  fullfile(xan_db_dir,'defaults','eeg');
        output3 =  fullfile(xan_db_dir,'defaults','meg');
    case 'xan_db_dir'
        output1 = fullfile(getUserDir(),'.XiAlphaNet','Datasets');
    case 'datasets'
        xan_db_dir = xan_get('xan_db_dir');        
        datasets_file =  fullfile(xan_db_dir,'Datasets.json');
        if(isfile(datasets_file))
            output1 = jsondecode(fileread(datasets_file));
        else
            output1 = [];
        end
    case 'dataset_info'
        datasets = xan_get('datasets');
        output1 = datasets(find(ismember({datasets.UUID},data),1));
        output2 = fullfile(output1.Location,'XIALPHANET.json');
    case 'datasets_file'
        output1 = fullfile(xan_get('xan_dir'),'Datasets','Datasets.json');
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
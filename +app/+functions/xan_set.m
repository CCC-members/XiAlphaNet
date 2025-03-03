function xan_set( varargin )
import app.*
import app.functions.*
import functions.*
import guide.*
import tools.*
if ((nargin >= 2) && ischar(varargin{1}))
    contextName = varargin{1};
    data = varargin{2};    
else
    return
end

% Get required context structure
switch contextName
    case 'xan_dir'
        
    case 'defaults_dir'
        
    case 'datasets'
        if(~isfolder(xan_get('xan_dir')))
            mkdir(xan_get('xan_dir'));
        end
        if(~isfolder(xan_get('xan_db_dir')))
            mkdir(xan_get('xan_db_dir'));
        end
        saveJSON(data,xan_get('datasets_file')) 
        h = matlab.desktop.editor.openDocument(xan_get('datasets_file'));
        h.smartIndentContents
        h.save
        h.close
    case 'dataset_info' 
        [~,file_name] = xan_get('dataset_info',data.UUID);
        saveJSON(data,file_name);
        h = matlab.desktop.editor.openDocument(file_name);
        h.smartIndentContents
        h.save
        h.close
    case 'tmp_path'
        properties = get_properties();
        properties.general_params.tmp.path = data;         
        
    case 'openmeeg'
        
    case 'duneruro'
end

end
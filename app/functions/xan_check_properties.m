function checked = xan_check_properties( varargin )
checked = true;
if ((nargin >= 1) && ischar(varargin{1}))
    contextName = varargin{1};
else
    return
end
if((nargin >= 2))
    path = varargin{2};
    path = strrep(path,'\','/');
else
    return;
end
if((nargin >= 3))
    ref_file = varargin{3};
end
% Check case parameters
switch contextName
    %% General
    case 'Input_path'
        if(~isfolder(fullfile(path)) && ~isempty(path))
            checked = false;
            return;
        end
    case 'Ref_file'
        if(~isfolder(path))
            checked = false;
            return;
        end
        structures = dir(path);
        structures(ismember( {structures.name}, {'.', '..'})) = [];  %remove . and ..
        count = 0;
        for i=1:length(structures)
            structure = structures(i);
            SubID = structure.name;
            tmp_file = strrep(ref_file,'SubID',SubID);
            if(~isfile(fullfile(structure.folder,structure.name,tmp_file)))
                count = count + 1;
            end
        end
        if(isequal(count,length(structures)))
            checked = false;
        end
    case 'Tmp_path'
        if(~isfolder(fullfile(path)) && ~isequal(path,'local'))
            checked = false;
            return;
        end
        if(isequal(path,'local'))
            return;
        end
        [~,values] = fileattrib(path);
        if(~values.UserWrite)
            checked = false;
        end
    
    %% Anatomy
    case 'Output_path'
        if(~isfolder(fullfile(path)) && ~isempty(path))
            checked = false;
            return;
        end
        [~,values] = fileattrib(path);
        if(~values.UserWrite)
            checked = false;
        end
    case 'Cortex'
        if(~isfile(path))
             checked = false;
             return;
        end
        
    case 'Leadfield'
        if(~isfile(path))
             checked = false;
             return;
        end
        
    case 'Channels'
        if(~isfile(path))
            checked = false;
            return;
        end        
  
    case 'Anat_conn'
        

    case 'Neuro_trac'

    case 'Cond_delay'

  

    
    otherwise
        checked = false;
end

end


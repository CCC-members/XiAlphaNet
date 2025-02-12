function [status, properties] = check_properties(properties)
%CHECK_PROPERTIES Summary of this function goes here
%   Detailed explanation goes here
status = true;
disp("-->> Checking properties");

if(~isfolder(properties.general_params.input_path))
    fprintf(2,strcat('\n-->> Error: The param input_path defined on app/general_params.json file: \n'));
    disp(properties.general_params.input_path);
    fprintf(2,strcat('It is not a correct adreess directory. \n'));
    disp('Please verify the location path.');
    status = false;
    return;
end
if(~isfolder(properties.general_params.output_path))
    fprintf(2,strcat('\n-->> Error: The param output_path defined on app/general_params.json file: \n'));
    disp(properties.general_params.output_path);
    fprintf(2,strcat('It is not a correct adreess directory. \n'));
    disp('Please verify the location path.');
    status = false;
    return;
else
    [status,values] = fileattrib(properties.general_params.output_path);
    if(~values.UserWrite)
        fprintf(2,strcat('\n-->> Error: The current user do not have write permissions on: \n'));
        disp(properties.general_params.output_path);
        disp('Please check the folder permission.');
        status = false;
        return;
    end
end

tmp_path = properties.general_params.tmp.path;
if(isequal(tmp_path,'local'))
    tmp_path = fullfile(pwd,'tmp');
    if(~isfolder(tmp_path))
        mkdir(tmp_path);
    end
    properties.general_params.tmp.path = tmp_path;
else
    if(~isfolder(tmp_path))
        fprintf(2,strcat('\n-->> Error: The param tmp.path defined on app/general_params.json file: \n'));
        disp(tmp_path);
        fprintf(2,strcat('It is not a correct adreess directory. \n'));
        disp('Please verify the location path.');
        status = false;
        return;
    else
        [status,values] = fileattrib(tmp_path);
        if(~values.UserWrite)
            fprintf(2,strcat('\n-->> Error: The current user do not have write permissions on: \n'));
            disp(tmp_path);
            disp('Please check the folder permission.');
            status = false;
            return;
        end
    end
end

anatomy = properties.anatomy_params;
if(~isfile(anatomy.leadfield.file_name))
    fprintf(2,strcat('\n-->> Error: The param leadfield.file_name defined on app/anatomy_params.json file: \n'));
    disp(anatomy.leadfield.file_name);
    fprintf(2,strcat('It is not a correct file. \n'));
    disp('Please verify the file location.');
    status = false;
    return;
end

if(~isfile(anatomy.cortex.file_name))
    fprintf(2,strcat('\n-->> Error: The param cortex.file_name defined on app/anatomy_params.json file: \n'));
    disp(anatomy.cortex.file_name);
    fprintf(2,strcat('It is not a correct file. \n'));
    disp('Please verify the file location.');
    status = false;
    return;
end

if(~isfile(anatomy.conn_anat.file_name))
    fprintf(2,strcat('\n-->> Error: The param conn_matrix.file_name defined on app/anatomy_params.json file: \n'));
    disp(anatomy.conn_matrix.file_name);
    fprintf(2,strcat('It is not a correct file. \n'));
    disp('Please verify the file location.');
    status = false;
    return;
end

if(~isfile(anatomy.cond_delay.file_name))
    fprintf(2,strcat('\n-->> Error: The param cond_delay.file_name defined on app/anatomy_params.json file: \n'));
    disp(anatomy.cond_delay.file_name);
    fprintf(2,strcat('It is not a correct file. \n'));
    disp('Please verify the file location.');
    status = false;
    return;
end

if(~isfile(anatomy.cond_delay.file_name))
    fprintf(2,strcat('\n-->> Error: The param cond_delay.file_name defined on app/anatomy_params.json file: \n'));
    disp(anatomy.cond_delay.file_name);
    fprintf(2,strcat('It is not a correct file. \n'));
    disp('Please verify the file location.');
    status = false;
    return;
end

if(~isfile(anatomy.channel.file_name))
    fprintf(2,strcat('\n-->> Error: The param channel.file_name defined on app/anatomy_params.json file: \n'));
    disp(anatomy.channel.file_name);
    fprintf(2,strcat('It is not a correct file. \n'));
    disp('Please verify the file location.');
    status = false;
    return;
end
end


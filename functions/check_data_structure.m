function [data,status,errors] = check_data_structure(properties, data)

Nw = properties.general_params.data.nFreqs;
errors = {};
status = true;

if(~isfield(data,'age'))
    data.age = randi([20,80],1);
end
if ischar(data.age) || isstring(data.age)
    % Convert the string to a number
    data.age = str2double(data.age);
elseif(isnan(data.age))
    data.age = randi([20 80],1);
else
    data.age = data.age;
end

if(~isfield(data,'freqrange'))
    status = false;
end

if(~isfield(data,'CrossM'))
    status = false;
else   
    % Cross
    data.Cross = data.CrossM(:,:,1:Nw);
    data.Cross = aveReference(data.Cross);
    data.freq = data.freqrange(1:Nw); 
    data = rmfield(data,{'CrossM','freqrange'});
end



end


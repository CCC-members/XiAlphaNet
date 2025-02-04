function saveJSON(data,output_file)
%SAVEJSON Summary of this function goes here
%   Detailed explanation goes here

data = jsonencode(data);
semi_idx = find(data == ',');
for i=1:length(semi_idx)
    if(isequal(data(semi_idx(i)+1),'"'))
        data(semi_idx(i)) = sprintf('*');
    end
end
data = strrep(data, '*', sprintf(',\r\t'));
data = strrep(data, '[', sprintf('[\r\t'));
data = strrep(data, '{', sprintf('{\r\t'));
data = strrep(data, '}', sprintf('\r}'));
data = strrep(data, ']', sprintf('\r]'));
fid = fopen(output_file, 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, data, 'char');
fclose(fid);
end

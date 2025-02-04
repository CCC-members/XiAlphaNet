function status = set_properties(properties)
%SET_PROPERTIES Summary of this function goes here
%   Detailed explanation goes here
status = true;

saveJSON(properties,strcat('app/properties.json'));


end


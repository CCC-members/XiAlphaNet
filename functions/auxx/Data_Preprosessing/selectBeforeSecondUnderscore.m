function result = selectBeforeSecondUnderscore(inputStr)
    % Find the positions of underscores in the input string
    underscorePositions = strfind(inputStr, '_');
    
    % Check if there are at least two underscores
    if length(underscorePositions) < 2
        error('Input string does not contain at least two underscores');
    end
    
    % Find the position of the second underscore
    secondUnderscorePos = underscorePositions(2);
    
    % Extract the substring before the second underscore
    result = inputStr(1:secondUnderscorePos-1);
end


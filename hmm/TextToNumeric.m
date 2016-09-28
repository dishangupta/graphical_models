% This function converts the data from text to 
% numeric form.

function [X] = TextToNumeric(filename)

    % Enumarate observation space
    symbols=['A' 'E' 'I' 'O' 'U' ' ' 'B' 'C' 'D' 'F' 'G' 'H' 'J' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'X' 'Y' 'Z'];

    fileID = fopen(filename, 'r');
    line = fgets(fileID);
    T = size(line, 2);
    O = size(symbols, 2);
    sequence = zeros(1, T);
    X = zeros(1, T);
    
    % Store character sequence
    for t = 1:size(line,2)
        sequence(t) = line(1, t);
    end
    
    % Convert characters to numbers
    for o = 1:O
        char = symbols(o);
        X(find(sequence == char)) = o;
    end
    
    X = X(1:end-1); % Remove final space
    X = X(find(X ~= 0)); % Remove out of bounds symbols
    fclose(fileID);
end


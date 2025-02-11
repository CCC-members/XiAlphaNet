function [T] = source_T(parameters)
    % Define the directory where you want to save the files
    saveDir = 'C:\Users\pedro\Desktop\Phy-SCE\Data\T_source';  % Adjust the path as needed

    % Check if the directory exists, if not, create it
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    Nr = parameters.Dimensions.Nr;
    D = sparse(parameters.Model.D);
    C = sparse(nulldiag(parameters.Model.C));

    freq = parameters.Data.freq;
    Nsw = length(freq);

    I2 = speye(Nr); % Sparse identity matrix
 

    for j = 1:Nsw
       fprintf('--> Source Operator %d\n',j)
        omega = freq(j);
        % Compute T_omega using sparse matrices
        T_omega = I2 - C .* exp(-2 * pi * 1i * omega * D);
        T = inv(T_omega); % Store result in a cell array

        % Save each computed matrix to the specified directory
        filename = fullfile(saveDir, sprintf('T_slice_%d.mat', j));
        save(filename, 'T');

        clear T; % Clear the variable to save memory
    end
end

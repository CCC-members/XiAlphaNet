function [RinvT] = source_RT(parameters)
    % Define the directory where you want to save the files
    saveDir = 'C:\Users\pedro\Desktop\Phy-SCE\Data\T_source';  % Adjust the path as needed

    % Check if the directory exists, if not, create it
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    N= 360;
    Nr = parameters.Dimensions.Nr;
    D = parameters.Model.D;
    %L =sparse(parameters.Model.L);
    C = sparse(nulldiag(parameters.Model.C));
    R = voxel_roi_map;
    freq = parameters.Data.freq;
    Nsw = length(freq);

    I2 = speye(Nr); % Sparse identity matrix
    RinvT = zeros(N,Nr,Nsw);

    for j = 1:Nsw
       fprintf('--> Source Operator %d\n',j)
        omega = freq(j);
        % Compute T_omega using sparse matrices
        T_omega = I2 - C .* exp(-2 * pi * 1i * omega * D);
        %T_omega = I2 - C .* exp(-2 * pi * 1i * omega * L/30000);
        RinvT(:,:,j) = R/(T_omega);
    end
end

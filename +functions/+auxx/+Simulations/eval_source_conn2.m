function [c] = eval_source_conn2(x, parameters)
    % Extract dimensions and model parameters
   % Ne = parameters.Dimensions.Ne;
   % Nr = parameters.Dimensions.Nr;
   % D = parameters.Model.D;
   % C = nulldiag(parameters.Model.C);
    %R = voxel_roi_map();
    freq = parameters.Data.freq;
    Nsw = length(freq);
   % Sw = parameters.Data.Cross;
    %T = parameters.Model.T;
    T = parameters.Model.T;
    %K=parameters.Model.K;
   % R = voxel_roi_map;
    % Unpack the input vector x into e, a, and sigma^2 components
    [e, a, ~] = x2v(x);
    %R = R*pinv(K);

    % Iterate over each frequency and compute cxi and calpha
    for j = 1:Nsw
        %fprintf('Processing frequency %d of %d\n', j, Nsw);
        omega = freq(j);

        % Calculate xi_omega and alpha_omega based on the model parameters
        xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
        alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);

        % Compute the transformed covariance matrices for xi and alpha
        c(:,:,j) = computeTDT(T(:,:,j), xi_omega+alpha_omega);
    end
end

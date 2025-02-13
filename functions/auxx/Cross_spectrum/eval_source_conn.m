function [solution] = eval_source_conn(x, freq,T,K,R,properties)
    % Extract dimensions and model parameters
    conn_delay = properties.general_params.parallel.conn_delay;
    [~,Nv,Nsw]=size(T);
    [Nr,~] = size(R);
    K_inv = pinv(K) ;
    U = R*K_inv;
    % Unpack the input vector x into e, a, and sigma^2 components
    [e, a, ~] = x2v(x);
    %R = R*pinv(K);
    c = zeros(Nr,Nr,Nsw);
    act = zeros(Nv,Nsw);
    % Iterate over each frequency and compute cxi and calpha
    tic
    if conn_delay == 1
        parfor j = 1:Nsw
            %fprintf('Processing frequency %d of %d\n', j, Nsw);
            omega = freq(j);
            Tj_cross = U*T(:,:,j);
            Tj_act = K_inv*T(:,:,j);
            % Calculate xi_omega and alpha_omega based on the model parameters
            xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
            alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
            sources_par = xi_omega+alpha_omega;
            % Compute the transformed covariance matrices for xi and alpha
            c(:,:,j) = (computeTDT(Tj_cross, sources_par));
            act(:,j) = Tj_act*sources_par;
        end
    else
        for j = 1:Nsw
            %fprintf('Processing frequency %d of %d\n', j, Nsw);
            omega = freq(j);
            Tj_cross = U*T(:,:,j);
            Tj_act = K_inv*T(:,:,j);
            % Calculate xi_omega and alpha_omega based on the model parameters
            xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
            alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
            sources_par = xi_omega+alpha_omega;
            % Compute the transformed covariance matrices for xi and alpha
            c(:,:,j) = (computeTDT(Tj_cross, sources_par));
            act(:,j) = Tj_act*sources_par;
        end
    end
    toc
    solution.Cross = c;
    solution.Activations = act;
end

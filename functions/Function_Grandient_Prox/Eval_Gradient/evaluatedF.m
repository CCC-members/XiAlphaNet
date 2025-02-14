function  [dF,F,smoothF] = evaluatedF(x,Ne,Nv,T,sw,sp,nsf_band,Sw)%parameters)  
    
    % Nw = parameters.Dimensions.Nw;
    % Ne = parameters.Dimensions.Ne;
    % Nv = parameters.Dimensions.Nv;
    % T = parameters.Model.T;
    % sw = parameters.Stochastic.Sampled.sw;    % stoc sampled value of freq
    % sp = parameters.Stochastic.Sampled.sp;           % stoch sampled position of freq
    % nsf_band = parameters.Stochastic.Nsfreq; % number of freq stoch sampled in a band
    Nsw = length(sw(1,:));
    % Sw = parameters.Data.Cross;
    
    % Unpack x into e, a, and sigma^2
    [e,a,sigma2] = x2v(x);
    
    % Initialize the third term
    term3 = 0;
    dFs2 = 0;
    % Define I, C, L, K (you'll need actual values for these)
    I = eye(Ne);
    I2 = eye(Nv);
    % Loop over omega (frequency components)
    dXi_de = zeros(Nv, 3, Nsw);
    dAlpha_da = zeros(Nv, 4, Nsw);
    TOISOT = zeros(Nv,Nsw);
    for j = 1:Nsw
        omega = sw(1,j); % sampled frequency 
        p  = sp(j); % sampled freq position

        % Evaluate spectral density matrix
        S_omega = Sw(:,:,p);

        % Compute T_omega
        T_omega = T(:,:,p);

        % Compute xi_omega and alpha_omega
        xi_omega = e(:,1) ./ (1 + e(:,2) .* omega^2).^e(:,3);
        alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);

        % Compute derivatives
        [dXi_de(:, :, j), dAlpha_da(:, :, j)] = dXiAlpha_deda(e, a, omega);

        % Construct Sigma_omega
        Sigma_omega = sigma2 * I +  computeTDT(T_omega, xi_omega + alpha_omega);
       
        % Regularize Sigma
        Sigma_omega = regularize_tensor(Sigma_omega);
        %
        SO_omega = S_omega * pinv(Sigma_omega);
        OISO = Sigma_omega \ (I-SO_omega);
        TOISOT(:,j) = computeDiagonalElements(T_omega',OISO);
        % Compute trace and determinant terms
        term3 = term3 + real(+logdet(Sigma_omega) + trace(SO_omega))* sw(2,j)/nsf_band;
        dFs2 = dFs2 +real(trace(OISO) )* sw(2,j)/nsf_band;
    end

    % Combine all terms
    F = term3;
    smoothF = term3;
    %% Gradient 
   dFa = zeros(Nv, 4);  % Gradient of the objective function with respect to the parameters of the alpha model
   dFe = zeros(Nv, 3); 

    for i = 1:Nv
        % Defining summation parameters for Alpha
        term1a = 0;
        term2a = 0;
        term3a = 0;
        term4a = 0;

        % Defining summation parameters for Xi
        term1e = 0;
        term2e = 0;
        term3e = 0; 

        for j = 1:Nsw
            %% Derivative with respect to vec
            df = real(TOISOT(i, j));

            %% Alpha
            % Alpha model derivative  
            edalpha_i1 = dAlpha_da(i, 1, j);
            edalpha_i2 = dAlpha_da(i, 2, j);
            edalpha_i3 = dAlpha_da(i, 3, j);
            edalpha_i4 = dAlpha_da(i, 4, j);

            % Summation over the frequencies
            term1a = term1a + df * edalpha_i1 * sw(2, j) / nsf_band;
            term2a = term2a + df * edalpha_i2 * sw(2, j) / nsf_band;
            term3a = term3a + df * edalpha_i3 * sw(2, j) / nsf_band;
            term4a = term4a + df * edalpha_i4 * sw(2, j) / nsf_band;

            %% Xi 
            edxi_i1 = dXi_de(i, 1, j); 
            edxi_i2 = dXi_de(i, 2, j);
            edxi_i3 = dXi_de(i, 3, j);

            % Summation over the frequencies
            term1e = term1e + df * edxi_i1 * sw(2, j) / nsf_band;
            term2e = term2e + df * edxi_i2 * sw(2, j) / nsf_band;
            term3e = term3e + df * edxi_i3 * sw(2, j) / nsf_band;
        end

        dFa(i, :) = [term1a, term2a, term3a, term4a];
        dFe(i, :) = [term1e, term2e, term3e];
    end
    % Vectorize gradients by row
    dFe_vec = dFe(:); % Vectorize dFe by row
    dFa_vec = dFa(:); % Vectorize dFa by row

    % Combine gradients
    dF = v2x(dFe_vec, dFa_vec, dFs2);
  end


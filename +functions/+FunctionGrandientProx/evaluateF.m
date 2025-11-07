function [F] = evaluateF(x,Ne,T,sw,sp,nsf_band,Sw)%;parameters)

import functions.*
import functions.auxx.*
import functions.auxx.ModelVectorization.*
import functions.auxx.OptimizedOperations.*

% Ne = parameters.Dimensions.Ne;
% T = parameters.Model.T;
% sw = parameters.Stochastic.Sampled.sw;    % stoc sampled value of freq
% sp = parameters.Stochastic.Sampled.sp;           % stoch sampled position of freq
% nsf_band = parameters.Stochastic.Nsfreq; % number of freq stoch sampled in a band
Nsw = length(sw(1,:));
% Sw = parameters.Data.Cross;
% Unpack x into e, a, and sigma^2
% Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[e,a,sigma2]  =x2v(x);

% Initialize the third term
term3 = 0;

% Define I
I = eye(Ne);

% Loop over omega (frequency components)
for j = 1:Nsw
    omega = sw(1,j); % sampled frequency
    p  = sp(j); % sampled freq position

    % Evaluate spectral density matrix
    S_omega = Sw(:,:,p);

    % Compute T_omega
    T_omega = T(:,:,p);

    % Compute xi_omega and alpha_omega
    xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
    alpha_omega = 0.5 * ( ...
        a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3) + ...
        a(:,1) ./ (1 + a(:,2) .* ((-omega) - a(:,4)).^2).^a(:,3) );
    % Define Sigma_omega
    Sigma_omega = sigma2 * I +  computeTDT(T_omega, xi_omega + alpha_omega);
    % Regularize Sigma
    Sigma_omega = regularize_tensor(Sigma_omega);

    % Compute trace and determinant terms
    term3 = term3 + real(+logdet(Sigma_omega) + trace(S_omega*pinv(Sigma_omega)))* sw(2,j)/nsf_band;
end


% Combine all terms
F = term3;
end

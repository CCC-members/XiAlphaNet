function [Sigma_omega] = reconstructSDM(x,parameters,freq)
import functions.auxx.StochasticEval.*
import functions.auxx.GenerateSourceSample.*
import functions.FunctionGrandientProx.*
import functions.auxx.ModelVectorization.*
import functions.auxx.OptimizedOperations.*
Ne = parameters.Dimensions.Ne;
T = parameters.Model.T;
parameters.Stochastic.stoch = 0;
[~,sw,sp] = sample_frequencies(freq,parameters.Stochastic.stoch,47);
Nsw = length(sw(1,:));
%Sw = data.Cross;
% Unpack x into e, a, and sigma^2
% Assume x is structured as [vec_r(e); vec_r(a); sigma^2]
[e,a,sigma2]  =x2v(x);

% Define I
I = eye(Ne);

% Loop over omega (frequency components)
for j = 1:Nsw
    omega = sw(1,j); % sampled frequency
    p  = sp(j); % sampled freq position

    % Evaluate spectral density matrix
    %S_omega = Sw(:,:,p);

    % Compute T_omega
    T_omega = T(:,:,p);

    % Compute xi_omega and alpha_omega
    xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
    alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);

    % Define Sigma_omega
    S = sigma2* I +  computeTDT(T_omega, xi_omega + alpha_omega);

    % Perform Singular Value Decomposition (SVD)
    % [U, S, V] = svd(S);

    % Regularize Sigma by ensuring the singular values are not too close to zero
    %tol = 1e-3; % Regularization tolerance
    %S = diag(max(real(diag(S)), tol));

    % Reconstruct the regularized Sigma
    Sigma_omega(:,:,j) =  S;

end

end
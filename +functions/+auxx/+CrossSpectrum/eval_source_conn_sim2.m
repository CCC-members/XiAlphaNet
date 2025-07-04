
function [Spectra] = eval_source_conn_sim2(x, freq,parameters)
import functions.*
import functions.auxx.*
import functions.auxx.Regularization.*
import functions.auxx.OptimizedOperations.*
import functions.auxx.ModelVectorization.*

% Extrat Dimensions
Nw = parameters.Dimensions.Nw;
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
%
if Nv==Nr
    C= nulldiag(parameters.Compact_Model.C);
    D = parameters.Compact_Model.D;
    R = eye(size(C));
else
    C = nulldiag(parameters.Model.C);
    D = parameters.Model.D;
end
C = 0.5*(C+C');
D = 0.5*(D+D');
I = speye(Nv);

% Unpack the input vector x into e, a, and sigma^2 components

[e, a, ~] = x2v(x);
spectra =zeros(Nv,Nw);
% Iterate over each frequency and compute cxi and calpha

for j = 1:Nw
    j
    %fprintf('Processing frequency %d of %d\n', j, Nsw);
    omega = freq(j);
    Z =  I +  C .* exp(-2 * pi * 1i * omega * D);
    %Tj_cross = R*Z;
    % Calculate xi_omega and alpha_omega based on the model parameters
    xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
    alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
    % Compute the transformed covariance matrices for xi and alpha
    spectra(:,j) = sum((abs(Z).^2) .* (xi_omega+alpha_omega), 2);
end

Spectra = spectra;
end


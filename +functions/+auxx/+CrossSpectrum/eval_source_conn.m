function [solution] = eval_source_conn(x, freq, R, properties, parameters)
import functions.*
import functions.auxx.*
import functions.auxx.Regularization.*
import functions.auxx.OptimizedOperations.*
import functions.auxx.ModelVectorization.*

% Extract dimensions and model parameters
conn_delay = properties.general_params.parallel.conn_delay;

% Extrat Dimensions
Nw = parameters.Dimensions.Nw;
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;

% Read Structural Data
if Nv==Nr
    c     = nulldiag(parameters.Compact_Model.C);
    D     = parameters.Compact_Model.D;
    R     = eye(size(C));
else
    C     = nulldiag(parameters.Model.C);
    D     = parameters.Model.D;
end
C         = 0.5*(C+C');
D         = 0.5*(D+D');
I         = speye(Nv);

% Unpack the input vector x into e, a, and sigma^2 components
[e, a, ~] = x2v(x);
c_full    = zeros(Nr,Nr,Nw);
c_xi      = zeros(Nr,Nr,Nw);
c_alpha   = zeros(Nr,Nr,Nw);
j_full    = zeros(Nv,Nw);
j_xi      = zeros(Nv,Nw);
j_alpha   = zeros(Nv,Nw);
s_full =zeros(Nv,Nw);

% Iterate over each frequency and compute cxi and calpha
if conn_delay == 1
    parfor j = 1:floor(Nw/2)
   
        %fprintf('Processing frequency %d of %d\n', j, Nsw);
        omega          = freq(j);
        Z              =  I +  C .* exp(-2 * pi * 1i * omega * D); % Approximate inverse
        Tj_cross       = R*Z;
        Tj_act         = Z;
        % Calculate xi_omega and alpha_omega based on the model parameters
        xi_omega       = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
        alpha_omega    = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
        %sources_par = xi_omega+alpha_omega;
        % Compute the transformed covariance matrices for xi and alpha
        c_xi(:,:,j)    = computeTDT(Tj_cross, xi_omega);
        c_alpha(:,:,j) = computeTDT(Tj_cross, alpha_omega);
        c_full(:,:,j)  = c_xi(:,:,j)+c_alpha(:,:,j) ;
        s_full(:,j)    = sum((abs(Z).^2) .* (xi_omega+alpha_omega), 2);
        j_xi(:,j)      = Tj_act*xi_omega;
        j_alpha(:,j)   = Tj_act*alpha_omega;
        j_full(:,j)    = j_xi(:,j)+j_alpha(:,j);
    end
     parfor j = floor(Nw/2)+1:Nw
        
        %fprintf('Processing frequency %d of %d\n', j, Nsw);
        omega          = freq(j);
        Z              =  I +  C .* exp(-2 * pi * 1i * omega * D); % Approximate inverse
        Tj_cross       = R*Z;
        Tj_act         = Z;
        % Calculate xi_omega and alpha_omega based on the model parameters
        xi_omega       = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
        alpha_omega    = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
        %sources_par = xi_omega+alpha_omega;
        % Compute the transformed covariance matrices for xi and alpha
        c_xi(:,:,j)    = computeTDT(Tj_cross, xi_omega);
        c_alpha(:,:,j) = computeTDT(Tj_cross, alpha_omega);
        c_full(:,:,j)  = c_xi(:,:,j)+c_alpha(:,:,j) ;
        s_full(:,j)    = sum((abs(Z).^2) .* (xi_omega+alpha_omega), 2);
        j_xi(:,j)      = Tj_act*xi_omega;
        j_alpha(:,j)   = Tj_act*alpha_omega;
        j_full(:,j)    = j_xi(:,j)+j_alpha(:,j);
    end
else
    for j = 1:Nw
        %fprintf('Processing frequency %d of %d\n', j, Nsw);
        omega          = freq(j);
        Z              =  I +  C .* exp(-2 * pi * 1i * omega * D);
        Tj_cross       = R*Z(:,:,j);
        Tj_act         = Z(:,:,j);
        % Calculate xi_omega and alpha_omega based on the model parameters
        xi_omega       = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
        alpha_omega    = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
        %sources_par = xi_omega+alpha_omega;
        % Compute the transformed covariance matrices for xi and alpha
        c_xi(:,:,j)    = computeTDT(Tj_cross, xi_omega);
        c_alpha(:,:,j) = computeTDT(Tj_cross, alpha_omega);
        c_full(:,:,j)  = c_xi(:,:,j)+c_alpha(:,:,j) ;
        s_full(:,j)    = sum((abs(Z).^2) .* (xi_omega+alpha_omega), 2);
        j_xi(:,j)      = Tj_act*xi_omega;
        j_alpha(:,j)   = Tj_act*alpha_omega;
        j_full(:,j)    = j_xi(:,j)+j_alpha(:,j);
    end
end

solution.Cross.Full        = c_full;
solution.Cross.Xi          = c_xi;
solution.Cross.Alpha       = c_alpha;
solution.Cross.Spectra     = s_full;
solution.Activations.Full  = j_full;
solution.Activations.Xi    = j_xi;
solution.Activations.Alpha = j_alpha;

end

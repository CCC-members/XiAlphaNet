function [solution] = eval_source_conn_sim(x, freq,G,K,R,properties)
import functions.*
import functions.auxx.*
import functions.auxx.OptimizedOperations.*
import functions.auxx.ModelVectorization.*

% Extract dimensions and model parameters
conn_delay = properties.general_params.parallel.conn_delay;
[~,Nv,Nsw]=size(G);
[Nr,~] = size(R);

% Unpack the input vector x into e, a, and sigma^2 components
[e, a, ~] = x2v(x);
%R = R*pinv(K);
c_full = zeros(Nr,Nr,Nsw);
c_xi = zeros(Nr,Nr,Nsw);
c_alpha = zeros(Nr,Nr,Nsw);
j_full = zeros(Nv,Nsw);
j_xi = zeros(Nv,Nsw);
j_alpha = zeros(Nv,Nsw);
% Iterate over each frequency and compute cxi and calpha
if conn_delay == 1
    parfor j = 1:Nsw
        %fprintf('Processing frequency %d of %d\n', j, Nsw);
        omega = freq(j);
        Tj_cross = G(:,:,j);
        Tj_act =   G(:,:,j);
        % Calculate xi_omega and alpha_omega based on the model parameters
        xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
        alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
        %sources_par = xi_omega+alpha_omega;
        % Compute the transformed covariance matrices for xi and alpha
        c_xi(:,:,j) = computeTDT(Tj_cross, xi_omega);
        c_alpha(:,:,j) = computeTDT(Tj_cross, alpha_omega);
        c_full(:,:,j) = c_xi(:,:,j)+c_alpha(:,:,j) ;
        j_xi(:,j) = Tj_act*xi_omega;
        j_alpha(:,j) = Tj_act*alpha_omega;
        j_full(:,j) = j_xi(:,j)+j_alpha(:,j);
    end
else
    for j = 1:Nsw
        %fprintf('Processing frequency %d of %d\n', j, Nsw);
        omega = freq(j);
        Tj_cross = G(:,:,j);
        Tj_act = G(:,:,j);
        % Calculate xi_omega and alpha_omega based on the model parameters
        xi_omega = e(:,1) ./ (1 + e(:,2) .* omega.^2).^e(:,3);
        alpha_omega = a(:,1) ./ (1 + a(:,2) .* (omega - a(:,4)).^2).^a(:,3);
        %sources_par = xi_omega+alpha_omega;
        % Compute the transformed covariance matrices for xi and alpha
        c_xi(:,:,j) = computeTDT(Tj_cross, xi_omega);
        c_alpha(:,:,j) = computeTDT(Tj_cross, alpha_omega);
        c_full(:,:,j) = c_xi(:,:,j)+c_alpha(:,:,j) ;
        j_xi(:,j) = Tj_act*xi_omega;
        j_alpha(:,j) = Tj_act*alpha_omega;
        j_full(:,j) = j_xi(:,j)+j_alpha(:,j);
    end
end

solution.Cross.Full = c_full;
solution.Cross.Xi = c_xi;
solution.Cross.Alpha = c_alpha;
solution.Activations.Full = j_full;
solution.Activations.Xi = j_xi;
solution.Activations.Alpha = j_alpha;
end

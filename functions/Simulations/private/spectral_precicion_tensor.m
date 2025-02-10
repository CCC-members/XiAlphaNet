function [sigma,theta,F] = spectral_precicion_tensor(theta0, nu0,  tstudent_a, tstudent_b, tstudent_dof, deltaf, Fmin, Fmax)
%% creationg a covariance/precision tensor with specific spectral form and assumptions of the underlying spectral factorization,
% the assumptions of the spectral factorization are encoded by the function "apply tensor phase"
%% Initializing
F           = Fmin:deltaf:Fmax;
q           = size(theta0,1);
Iq          = eye(q);
theta       = zeros(q,q,length(F));
sigma       = zeros(q,q,length(F));
mask0       = ones(q,q);
% mask0(abs(theta0) == 0) = 0;
%% crating tensors
disp(['<spectral_precicion_tensor> Generating ',num2str(length(F)),' slices cross spectrum;'])
phi = tstudent_spectra(F, tstudent_a,tstudent_b,tstudent_dof,nu0);
phi = phi/max(phi);%scaling to 0-1;
% phi(phi < 0.001) = 0.001;
K = apply_tensor_phase(theta0, nu0, deltaf, Fmin, Fmax);
K = K./repmat(reshape(phi,1,1,length(F)),q,q,1);
U = Iq - K;
for count = 1:length(F)
    if (7 <= F(count)) && (F(count) <= 14)
        theta(:,:,count) = (U(:,:,count)'*U(:,:,count)).*mask0;
        sigma(:,:,count) = inv(theta(:,:,count));
    else
        theta(:,:,count) = U(:,:,count)'*U(:,:,count);
        sigma(:,:,count) = inv(theta(:,:,count));
    end
end
%% set DC coef as zero
theta(:,:,1)=0;
sigma(:,:,1)=0;
end

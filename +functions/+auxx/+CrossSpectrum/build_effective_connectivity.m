function [Cxi, Calpha] = build_effective_connectivity(subj_folder, params_folder,R)
% BUILD_EFFECTIVE_CONNECTIVITY
%   Construct Xi and Alpha effective connectivity for one subject
%
% INPUTS:
%   subj_folder   - path to subject’s folder (with Alpha_estimate.mat, Xi_estimate.mat, Mod_Weights.mat)
%   params_folder - path to folder containing parameters.mat
%   freq_idx_xi   - index of frequency bin for Xi process
%   freq_idx_alpha- index of frequency bin for Alpha process
%
% OUTPUTS:
%   Cxi     - Xi effective connectivity matrix
%   Calpha  - Alpha effective connectivity matrix
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*
import functions.auxx.OptimizedOperations.*
% --- Load subject spectral parameters ---
A = load(fullfile(subj_folder, "Alpha_estimate.mat"));
X = load(fullfile(subj_folder, "Xi_estimate.mat"));
MW = load(fullfile(subj_folder, "Mod_Weights.mat"));

% --- Load structural model ---
parameters = load(fullfile(params_folder, "parameters.mat"));
C0 = parameters.Model.C;
D0 = parameters.Model.D;
Nr = size(R,1);

% --- Apply subject-specific scaling ---
w_delay = MW.Mod_Weights(1);   % scale for delays
w_conn  = MW.Mod_Weights(2);   % scale for connectivity
C = 1  * C0;
D = 1 * D0;

omega_xi    = 0.3;
omega_alpha = 10;

% --- Build transfer operators ---
I = eye(size(C));
Z_xi    = I + C .* exp(-2*pi*1i*omega_xi*D);
Z_alpha = I + C .* exp(-2*pi*1i*omega_alpha*D);

T_xi    = R * Z_xi;
T_alpha = R * Z_alpha;

% --- Spectral weights ---
xi_omega = X.Power ./ (1 + X.Width .* omega_xi.^2).^X.Exponent;

alpha_omega = 0.5 * ( ...
    A.Power ./ (1 + A.Width .* (omega_alpha - A.PAF).^2).^A.Exponent + ...
    A.Power ./ (1 + A.Width .* (omega_alpha + A.PAF).^2).^A.Exponent );

% --- Effective connectivity matrices ---
Cxi    = computeTDT(T_xi, xi_omega);
Calpha = computeTDT(T_alpha, alpha_omega);

end

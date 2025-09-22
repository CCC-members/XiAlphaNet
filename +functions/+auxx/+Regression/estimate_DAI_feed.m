clc;
clear all;

% === Imports ===
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*
import functions.auxx.OptimizedOperations.*

% === Parameters ===
json_path = 'E:\test262025_delete\XIALPHANET.json';
load("+templates\Cortex_with_myelin.mat")
% hierarchy.mat must contain variable h (Nr x 1 vector of hierarchy scores)

% Load hierarchy scores
h = Cortex.Atlas(Cortex.iAtlas).MyelinValues;   % one hierarchy value per ROI

% === Load dataset metadata ===
[dataset_dir, ~, ~] = fileparts(json_path);
dataset = jsondecode(fileread(json_path));
dataset.Location = dataset_dir;

parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));
R = parameters.Model.R;
C0 = parameters.Model.C;
D0 = parameters.Model.D;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
Nw = parameters.Dimensions.Nw;

% === Frequency indices for Xi and Alpha ===
% Set data dir
dir_data = 'E:\norms';
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
selected_folders = subject_folders(randperm(length(subject_folders), 1));
subject_folder = selected_folders(1).name;
mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
data_struct = load(mat_file_path);
freq = data_struct.data_struct.freqrange(1:Nw);
parameters.Data.freq = freq;

%
omega_xi = parameters.Data.freq(1);    % low-frequency Xi
omega_alpha = parameters.Data.freq(25);% alpha peak example

% === Collect participants ===
n = length(dataset.Participants);
All_Data_tmp = cell(1, n);
Mod_Matrix = cell(1, n);

ages_tmp = nan(1, n);
for i = 1:n
    n-i
    participant = dataset.Participants(i);
    if isequal(participant.Status, 'Completed')
        try
            Part_Info = jsondecode(fileread(fullfile(dataset.Location, ...
                                   participant.SubID, participant.FileInfo)));

            Mod_Matrix{i} = load(fullfile(dataset.Location, ...
                                   participant.SubID, Part_Info.Mod_Weights));

            % Load Alpha
            alpha_process = load(fullfile(dataset.Location, ...
                                    participant.SubID, Part_Info.Alpha_estimate));
            a = [alpha_process.Power, ...
                 alpha_process.Width, ...
                 alpha_process.Exponent, ...
                 alpha_process.PAF];

            % Load Xi
            xi_process = load(fullfile(dataset.Location, ...
                                    participant.SubID, Part_Info.Xi_estimate));
            e = [xi_process.Power, ...
                 xi_process.Width, ...
                 xi_process.Exponent];

            All_Data_tmp{i} = {e,a};
            ages_tmp(i) = participant.Age;
        catch
            warning("Participant %s failed.", participant.SubID);
        end
    end
end

valid = ~cellfun(@isempty, All_Data_tmp);
All_Data = All_Data_tmp(valid);
Mod_Matrix = Mod_Matrix(valid);
ages = ages_tmp(valid);
N = length(All_Data);

% === Initialize accumulators ===
DAI_xi_all = zeros(Nr,Nr,N);
DAI_alpha_all = zeros(Nr,Nr,N);

% === Main loop ===
for s = 1:N
    s
    ea = All_Data{s};
    e = ea{1}; a = ea{2};

    % Subject-specific scaling
    C = Mod_Matrix{s}.Mod_Weights(2) * C0;
    D = Mod_Matrix{s}.Mod_Weights(1) * D0;
    I = eye(Nv);

    % Transfer operators
    Z_xi    = I + C .* exp(-2*pi*1i*omega_xi*D);
    Z_alpha = I + C .* exp(-2*pi*1i*omega_alpha*D);

    T_xi    = R * Z_xi;
    T_alpha = R * Z_alpha;

    % Spectral weights
    xi_omega = e(:,1) ./ (1 + e(:,2).*omega_xi.^2).^e(:,3);
    alpha_omega = 0.5*( a(:,1)./(1 + a(:,2).*(omega_alpha - a(:,4)).^2).^a(:,3) + ...
                        a(:,1)./(1 + a(:,2).*(omega_alpha + a(:,4)).^2).^a(:,3) );

    % Cross-spectral density matrices
    A_xi    = computeTDT(T_xi, xi_omega);
    A_alpha = computeTDT(T_alpha, alpha_omega);

    % Precision matrices
    P_xi    = pinv(A_xi);
    P_alpha = pinv(A_alpha);

    % Granger causality matrices
    GC_xi    = zeros(Nr);
    GC_alpha = zeros(Nr);
    for i = 1:Nr
        for j = 1:Nr
            if i==j, continue; end
            % Xi process
            pii = real(P_xi(i,i)); pjj = real(P_xi(j,j)); pij = real(P_xi(i,j));
            Sii_full = 1/pii;
            Sii_cond = 1/(pii - (pij^2)/pjj);
            if Sii_full>0 && Sii_cond>0
                GC_xi(j,i) = log(Sii_full/Sii_cond);
            end
            % Alpha process
            pii = real(P_alpha(i,i)); pjj = real(P_alpha(j,j)); pij = real(P_alpha(i,j));
            Sii_full = 1/pii;
            Sii_cond = 1/(pii - (pij^2)/pjj);
            if Sii_full>0 && Sii_cond>0
                GC_alpha(j,i) = log(Sii_full/Sii_cond);
            end
        end
    end

    % Directed Asymmetry Index (DAI)
    DAI_xi = zeros(Nr);
    DAI_alpha = zeros(Nr);
    for i = 1:Nr
        for j = 1:Nr
            if i==j, continue; end
            num_xi = GC_xi(i,j) - GC_xi(j,i);
            den_xi = GC_xi(i,j) + GC_xi(j,i) + eps;
            DAI_xi(i,j) = sign(h(i)-h(j)) * (num_xi/den_xi);

            num_alpha = GC_alpha(i,j) - GC_alpha(j,i);
            den_alpha = GC_alpha(i,j) + GC_alpha(j,i) + eps;
            DAI_alpha(i,j) = sign(h(i)-h(j)) * (num_alpha/den_alpha);
        end
    end

    DAI_xi_all(:,:,s) = DAI_xi;
    DAI_alpha_all(:,:,s) = DAI_alpha;
end

% === Population averages ===
DAI_xi_pop = mean(DAI_xi_all,3,'omitnan');
DAI_alpha_pop = mean(DAI_alpha_all,3,'omitnan');



%%  Hierh_DAI2

clc;
clear all;

% === Imports ===
import functions.auxx.ModelVectorization.*
import guide.Visualization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.Refine_Solution.*
import functions.auxx.OptimizedOperations.*

% Path to the JSON file with model result metadata. Modify this directions
% manually acording to the location of the downloaded data 
json_path = '/mnt/Develop/Ronaldo/program_working/xialphanet_newresults22/XIALPHANET.json';
dir_data = '/mnt/Develop/Ronaldo/dev/Data/norms';

% === Parameters ===
load("+templates/Cortex_with_myelin.mat")
% hierarchy.mat must contain variable h (Nr x 1 vector of hierarchy scores)

% Load hierarchy scores
h =zscore(Cortex.Atlas(Cortex.iAtlas).MyelinValues);   % one hierarchy value per ROI

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
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
selected_folders = subject_folders(randperm(length(subject_folders), 1));
subject_folder = selected_folders(1).name;
mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
data_struct = load(mat_file_path);
freq = data_struct.data_struct.freqrange(1:Nw);
parameters.Data.freq = freq;

omega_xi = parameters.Data.freq(1);      % low-frequency Xi
omega_alpha = parameters.Data.freq(25);  % alpha peak example

% === Collect participants ===
age_min = 0;
age_max = 100;

n = length(dataset.Participants);
All_Data_tmp = cell(1, n);
Mod_Matrix = cell(1, n);
ages_tmp = nan(1, n);

parfor i = 1:n
    i
    participant = dataset.Participants(i);
    participant_age = participant.Age;

    if isequal(participant.Status, 'Completed') && ...
            age_min <= participant_age && participant_age <= age_max
        try
            % Load participant info
            Part_Info = jsondecode(fileread(fullfile(dataset.Location, ...
                participant.SubID, participant.FileInfo)));

            % Load Mod Matrix
            Mod_Matrix{i} = load(fullfile(dataset.Location, ...
                participant.SubID, Part_Info.Mod_Weights));

            % === Alpha ===
            alpha_process = load(fullfile(dataset.Location, ...
                participant.SubID, Part_Info.Alpha_estimate));
            a = [alpha_process.Power, ...
                alpha_process.Width, ...
                alpha_process.Exponent, ...
                alpha_process.PAF];

            % Threshold on Power
            ta = set_threshold_em(a(:,1));
            a(:,1) = a(:,1) .* (a(:,1) > ta);

            % === Xi ===
            xi_process = load(fullfile(dataset.Location, ...
                participant.SubID, Part_Info.Xi_estimate));
            e = [xi_process.Power, ...
                xi_process.Width, ...
                xi_process.Exponent];

            % Threshold on Power
            te = set_threshold_em(e(:,1));
            e(:,1) = e(:,1) .* (e(:,1) > te);

            % Assign combined
            All_Data_tmp{i} = {e, a};
            ages_tmp(i) = participant_age;

        catch ME
            warning("Participant %s failed: %s", participant.SubID, ME.message);
        end
    end
end

% === Filter valid entries ===
valid = ~cellfun(@isempty, All_Data_tmp);
All_Data = All_Data_tmp(valid);
Mod_Matrix = Mod_Matrix(valid);
ages = ages_tmp(valid);
N = length(All_Data);
%% === Geweke spectral GC directly from Xi–AlphaNET transfer ===
fprintf('--- Stage 1: Averaging sqrt_xi, sqrt_alpha, and Mod_Weights ---\n');
N = numel(All_Data);
sqrt_xi_sum    = zeros(Nv,1);
sqrt_alpha_sum = zeros(Nv,1);
wD_sum = 0;
wC_sum = 0;
count_valid = 0;
batch_size  = 20;
num_batches = ceil(N / batch_size);
for b = 1:num_batches
    fprintf('Batch %d/%d (averaging spectral and structural params)...\n', b, num_batches);

    idx_start = (b-1)*batch_size + 1;
    idx_end   = min(b*batch_size, N);
    idx_batch = idx_start:idx_end;

    for k = 1:numel(idx_batch)
        s = idx_batch(k);
        ea = All_Data{s};
        e = ea{1}; a = ea{2};

        sqrt_xi    = sqrt(e(:,1) ./ (1 + e(:,2).*omega_xi.^2).^e(:,3));
        sqrt_alpha = sqrt(0.5*( a(:,1)./(1 + a(:,2).*(omega_alpha - a(:,4)).^2).^a(:,3) + ...
                                 a(:,1)./(1 + a(:,2).*(omega_alpha + a(:,4)).^2).^a(:,3) ));

        sqrt_xi_sum    = sqrt_xi_sum    + sqrt_xi;
        sqrt_alpha_sum = sqrt_alpha_sum + sqrt_alpha;

        wD_sum = wD_sum + Mod_Matrix{s}.Mod_Weights(1);
        wC_sum = wC_sum + Mod_Matrix{s}.Mod_Weights(2);

        count_valid = count_valid + 1;
    end
end

sqrt_xi_mean    = sqrt_xi_sum    / count_valid;
sqrt_alpha_mean = sqrt_alpha_sum / count_valid;
mean_wD = wD_sum / count_valid;
mean_wC = wC_sum / count_valid;

fprintf('--- Stage 2: Computing single averaged transfer functions ---\n');

I = eye(Nv);
C_mean = mean_wC * C0;
D_mean = mean_wD * D0;

inv_xi    = inv(I - C_mean .* exp(-2*pi*1i*omega_xi*D_mean));
inv_alpha = inv(I - C_mean .* exp(-2*pi*1i*omega_alpha*D_mean));

Hxi    = R*inv_xi    * diag(sqrt_xi_mean)*R';
Halpha = R*inv_alpha * diag(sqrt_alpha_mean)*R';

fprintf('--- Stage 3: Computing Geweke GC from averaged transfer functions ---\n');

Gxi_mean    = zeros(Nr);
Galpha_mean = zeros(Nr);

for i = 1:Nr
    Sii_xi    = sum(abs(Hxi(i,:)).^2);
    Sii_alpha = sum(abs(Halpha(i,:)).^2);

    for j = 1:Nr
        if i == j, continue; end
        Hij_xi    = Hxi(i,j);
        Hij_alpha = Halpha(i,j);

        denom_xi    = Sii_xi    - abs(Hij_xi)^2;
        denom_alpha = Sii_alpha - abs(Hij_alpha)^2;

        if denom_xi > 0
            Gxi_mean(i,j) = log(Sii_xi / denom_xi);
        end
        if denom_alpha > 0
            Galpha_mean(i,j) = log(Sii_alpha / denom_alpha);
        end
    end
end


%% === INPUT: Precomputed Granger causality matrices ===
G_xi    = Gxi_mean;
G_alpha = Galpha_mean;
Nr = size(G_xi, 1);
nPerm = 1000;
alpha_level = 0.01;

fprintf('Running %d permutations on GC for significance...\n', nPerm);

% === Step 1: Permutation null on GC ===
% Null hypothesis: directed flow pattern is random w.r.t. ROI indices
Gxi_null    = zeros(Nr, Nr, nPerm);
Galpha_null = zeros(Nr, Nr, nPerm);

for p = 1:nPerm
    perm = randperm(Nr);
    Gxi_null(:,:,p)    = G_xi(perm, perm);
    Galpha_null(:,:,p) = G_alpha(perm, perm);
end

% --- Compute empirical one-tailed p-values (large GC)
pG_xi    = zeros(Nr);
pG_alpha = zeros(Nr);

for i = 1:Nr
    for j = 1:Nr
        if i == j, continue; end
        null_xi = squeeze(Gxi_null(i,j,:));
        null_a  = squeeze(Galpha_null(i,j,:));
        pG_xi(i,j)    = mean(null_xi >= G_xi(i,j));
        pG_alpha(i,j) = mean(null_a  >= G_alpha(i,j));
    end
end

% --- Significance masks ---
mask_Axi    = pG_xi    < alpha_level;
mask_Aalpha = pG_alpha < alpha_level;

fprintf('Significant GC pairs: \\xi %.2f%%, \\alpha %.2f%%\n', ...
    100*mean(mask_Axi(:)), 100*mean(mask_Aalpha(:)) );

% === Step 2: Amplitude thresholding (optional) ===
thr_xi    = prctile(G_xi(:), 0.0);
thr_alpha = prctile(G_alpha(:), 0.0);

mask_amp_xi    = G_xi    > thr_xi;
mask_amp_alpha = G_alpha > thr_alpha;

fprintf('Amplitude thresholds: \\xi %.4f, \\alpha %.4f\n', thr_xi, thr_alpha);

mask_final_xi    = triu(mask_Axi    & mask_amp_xi, 1);
mask_final_alpha = triu(mask_Aalpha & mask_amp_alpha, 1);

[src_xi, tgt_xi]       = find(mask_final_xi);
[src_alpha, tgt_alpha] = find(mask_final_alpha);

fprintf('Selected %d \\xi-band and %d \\alpha -band GC connections.\n', ...
    numel(src_xi), numel(src_alpha));

%=== Step 3: Normalize GC values ===
GC_xi_norm    = G_xi / max(G_xi(:));
GC_alpha_norm = G_alpha / max(G_alpha(:));

fprintf('GC normalization completed.\n');

% === Step 4: Compute DAI restricted to significant GC ===
fprintf('Computing DAI restricted to significant GC...\n');
eps_val = 1e-30;
DAI_xi    = zeros(Nr);
DAI_alpha = zeros(Nr);

for i = 1:Nr
    for j = 1:Nr
        if i == j, continue; end

            num_xi = G_xi(i,j) - G_xi(j,i);
            den_xi = G_xi(i,j) + G_xi(j,i) + eps_val;
            DAI_xi(i,j) = sign(h(i) - h(j)) * num_xi / den_xi;

            num_a = G_alpha(i,j) - G_alpha(j,i);
            den_a = G_alpha(i,j) + G_alpha(j,i) + eps_val;
            DAI_alpha(i,j) = sign(h(i) - h(j)) * num_a / den_a;
    end
end

DAI_xi(~mask_final_xi)    = 0;
DAI_alpha(~mask_final_alpha) = 0;

DAI_xi    = DAI_xi / max(abs(DAI_xi(:)));
DAI_alpha = DAI_alpha / max(abs(DAI_alpha(:)));

fprintf('DAI computation completed.\n');



fprintf('DAI computed for significant GC connections only.\n');
 
% === Step 4. Visualization ===
figure('Color','w','Name','DAI from significant GC');
subplot(1,2,1);
imagesc(DAI_xi); axis square;
title('\xi-band DAI (significant connections only)','FontWeight','bold');
colorbar; xlabel('Source ROI'); ylabel('Target ROI');
 
subplot(1,2,2);
imagesc(DAI_alpha); axis square;
title('\alpha-band DAI (significant connections only)','FontWeight','bold');
colorbar; xlabel('Source ROI'); ylabel('Target ROI');
fprintf('--- DAI significance filtering complete ---\n');

target_DAI = DAI_alpha;
rotSteps = 0;
circos_plot_test;

target_DAI = DAI_xi;
rotSteps = -7;
circos_plot_test;

A = Gxi_mean;
guide.Visualization.effective_spectral_net;

A = Galpha_mean;
guide.Visualization.effective_spectral_net


%% === PARAMETERS ===
eps_val = 1e-12;

Nr = size(Gxi_mean,1);
DAI_xi    = zeros(Nr);
DAI_alpha = zeros(Nr);

% === Compute Directional Asymmetry Index (Only order defined by Myelin) ===
for i = 1:Nr
    for j = i+1:Nr
        % \xi-band
            val = sign(h(i) - h(j)) ;
            DAI_xi(i,j) = val;
            DAI_xi(j,i) = -val;

            val = sign(h(i) - h(j)) ;
            DAI_alpha(i,j) = val;
            DAI_alpha(j,i) = -val;
    end
end

% === Separate feedforward and feedback connections ===
feed_mask_xi     = DAI_xi > 0;
fb_mask_xi       = DAI_xi < 0;
feed_mask_alpha  = DAI_alpha > 0;
fb_mask_alpha    = DAI_alpha < 0;

feed_values_xi    = Gxi_mean(feed_mask_xi);
fb_values_xi      = Gxi_mean(fb_mask_xi);
feed_values_alpha = Galpha_mean(feed_mask_alpha);
fb_values_alpha   = Galpha_mean(fb_mask_alpha);

% === Observed statistics ===
feed_mean_xi  = mean(feed_values_xi,'omitnan');
fb_mean_xi    = mean(fb_values_xi,'omitnan');
diff_obs_xi   = fb_mean_xi - feed_mean_xi;

feed_mean_alpha = mean(feed_values_alpha,'omitnan');
fb_mean_alpha   = mean(fb_values_alpha,'omitnan');
diff_obs_alpha  = fb_mean_alpha - feed_mean_alpha;
% === Observed statistics ===
fprintf('=== Observed \\delta GC (Feedback - Feedforward) ===\n');
fprintf('\\xi-band:  \\delta GC = %.4f\n', diff_obs_xi);
fprintf('\\alpha-band:  \\delta GC = %.4f\n\n', diff_obs_alpha);

% === Build null distributions with balanced subsampling ===
nPerm = 1000;
subsample_n = 500;  % number of edges per group used for each permutation

all_vals_xi    = [feed_values_xi(:); fb_values_xi(:)];
n_feed_xi_full = numel(feed_values_xi);
n_fb_xi_full   = numel(fb_values_xi);
null_diff_xi   = zeros(nPerm,1);

all_vals_alpha    = [feed_values_alpha(:); fb_values_alpha(:)];
n_feed_alpha_full = numel(feed_values_alpha);
n_fb_alpha_full   = numel(fb_values_alpha);
null_diff_alpha   = zeros(nPerm,1);

for b = 1:nPerm
    % === ?-band ===
    % draw balanced random subset for permutation
    perm_idx = randperm(n_feed_xi_full + n_fb_xi_full);
    perm_vals = all_vals_xi(perm_idx);
    n_sub = min([subsample_n, floor(numel(perm_vals)/2)]);
    % split balanced subset into pseudo feed and feedback
    perm_feed = randsample(perm_vals, n_sub, false);
    perm_fb   = randsample(perm_vals, n_sub, false);
    null_diff_xi(b) = mean(perm_fb,'omitnan') - mean(perm_feed,'omitnan');

    % === a-band ===
    perm_idx = randperm(n_feed_alpha_full + n_fb_alpha_full);
    perm_vals = all_vals_alpha(perm_idx);
    n_sub = min([subsample_n, floor(numel(perm_vals)/2)]);
    perm_feed = randsample(perm_vals, n_sub, false);
    perm_fb   = randsample(perm_vals, n_sub, false);
    null_diff_alpha(b) = mean(perm_fb,'omitnan') - mean(perm_feed,'omitnan');
end

% === One-tailed permutation tests ===
p_xi_FBgtFF    = mean(null_diff_xi   >= diff_obs_xi);
p_xi_FFgtFB    = mean(null_diff_xi   <= diff_obs_xi);
p_alpha_FBgtFF = mean(null_diff_alpha>= diff_obs_alpha);
p_alpha_FFgtFB = mean(null_diff_alpha<= diff_obs_alpha);

fprintf('=== Permutation Tests ===\n');
fprintf('\\xi-band:  \\delta GC = %.4f  p_FB>FF = %.4f  p_FF>FB = %.4f\n', ...
        diff_obs_xi, p_xi_FBgtFF, p_xi_FFgtFB);
fprintf('\\alpha-band:  \\delta GC = %.4f  p_FB>FF = %.4f  p_FF>FB = %.4f\n\n', ...
        diff_obs_alpha, p_alpha_FBgtFF, p_alpha_FFgtFB);

% === Plotting (exact same style as reference) ===
figure('Color','w','Position',[100 100 1200 800],'Renderer','painters');
cmap = hot;  % black ? red ? yellow ? white

% ---- (1) ?-band: Feedback > Feedforward ----
subplot(2,2,1); hold on;
hh(1) = histogram(null_diff_xi,40,'Normalization','pdf','EdgeColor','none');
hh(1).FaceColor = cmap(round(end*0.8),:);
mu = mean(null_diff_xi); sigma = std(null_diff_xi);
xvals = linspace(min(null_diff_xi),max(null_diff_xi),200);
plot(xvals, normpdf(xvals,mu,sigma),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
xline(diff_obs_xi,'r','LineWidth',2);
title(['\xi-band: Feedback > Feedforward', newline, ...
       sprintf('\\DeltaGC = %.4f,  p = %.4f', diff_obs_xi, p_xi_FBgtFF)], ...
       'FontWeight','bold','FontSize',12);
xlabel('\DeltaGC = Feedback - Feedforward', 'FontSize',14,'FontWeight','bold');
ylabel('Density', 'FontSize',14,'FontWeight','bold');
ax_main = gca;
set(ax_main,'FontSize',13,'FontWeight','bold');
xlim([min(xvals) max(xvals)]); grid on; box off; axis tight;

% ---- (2) ?-band: Feedforward > Feedback ----
subplot(2,2,2); hold on;
hh(2) = histogram(null_diff_xi,40,'Normalization','pdf','EdgeColor','none');
hh(2).FaceColor = cmap(round(end*0.6),:);
plot(xvals, normpdf(xvals,mu,sigma),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
xline(diff_obs_xi,'r','LineWidth',2);
title(['\xi-band: Feedforward > Feedback', newline, ...
       sprintf('\\DeltaGC = %.4f,  p = %.4f', diff_obs_xi, p_xi_FFgtFB)], ...
       'FontWeight','bold','FontSize',12);
xlabel('\DeltaGC = Feedback - Feedforward', 'FontSize',14,'FontWeight','bold');
ylabel('Density', 'FontSize',14,'FontWeight','bold');
ax_main = gca;
set(ax_main,'FontSize',13,'FontWeight','bold');
xlim([min(xvals) max(xvals)]); grid on; box off; axis tight;

% ---- (3) a-band: Feedback > Feedforward ----
subplot(2,2,3); hold on;
hh(3) = histogram(null_diff_alpha,40,'Normalization','pdf','EdgeColor','none');
hh(3).FaceColor = cmap(round(end*0.8),:);
mu = mean(null_diff_alpha); sigma = std(null_diff_alpha);
xvals = linspace(min(null_diff_alpha),max(null_diff_alpha),200);
plot(xvals, normpdf(xvals,mu,sigma),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
xline(diff_obs_alpha,'r','LineWidth',2);
title(['\alpha-band: Feedback > Feedforward', newline, ...
       sprintf('\\DeltaGC = %.4f,  p = %.4f', diff_obs_alpha, p_alpha_FBgtFF)], ...
       'FontWeight','bold','FontSize',12);
xlabel('\DeltaGC = Feedback - Feedforward', 'FontSize',14,'FontWeight','bold');
ylabel('Density', 'FontSize',14,'FontWeight','bold');
ax_main = gca;
set(ax_main,'FontSize',13,'FontWeight','bold');
xlim([min(xvals) max(xvals)]); grid on; box off; axis tight;

% ---- (4) a-band: Feedforward > Feedback ----
subplot(2,2,4); hold on;
hh(4) = histogram(null_diff_alpha,40,'Normalization','pdf','EdgeColor','none');
hh(4).FaceColor = cmap(round(end*0.6),:);
plot(xvals, normpdf(xvals,mu,sigma),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
xline(diff_obs_alpha,'r','LineWidth',2);
title(['\alpha-band: Feedforward > Feedback', newline, ...
       sprintf('\\DeltaGC = %.4f,  p = %.4f', diff_obs_alpha, p_alpha_FFgtFB)], ...
       'FontWeight','bold','FontSize',12);
xlabel('\DeltaGC = Feedback - Feedforward', 'FontSize',14,'FontWeight','bold');
ylabel('Density', 'FontSize',14,'FontWeight','bold');
ax_main = gca;
set(ax_main,'FontSize',13,'FontWeight','bold');
xlim([min(xvals) max(xvals)]); grid on; box off; axis tight;

sgtitle('Permutation Tests for Directional Asymmetry (\DeltaGC)', ...
        'FontWeight','bold','FontSize',16);

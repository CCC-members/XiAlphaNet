clc; clear; close all;

%% CONFIG
root_dir       = "/mnt/Develop/Ronaldo/organized_data_test_retest/processed_data";
groups         = ["session1","session2"];
B_boot         = 1000;   % bootstrap resamples
B_perm         = 5000;   % permutations
rho0           = 0.4;    % null hypothesis threshold
q_FDR          = 0.05;   % FDR significance level
rng(0);

disp("-->> Starting process (ICC test + FDR)");

%% Load template cortex and ROI projection matrix
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit

Cortex   = load("templates/Cortex.mat");
template = load("templates/axes.mat");
currentAxes = template.axes;
hemis = [];
[Cortex, iHideVert] = split_hemisphere(Cortex, hemis);

% ROI projection matrix (Nv × 360)
[R,R_inv] = functions.auxx.ModelVectorization.roi_average_operators(Cortex,10);

%% Initialize subject maps
alpha_pre_map = containers.Map('KeyType','char','ValueType','any');
alpha_post_map= containers.Map('KeyType','char','ValueType','any');
xi_pre_map    = containers.Map('KeyType','char','ValueType','any');
xi_post_map   = containers.Map('KeyType','char','ValueType','any');

%% Load data
for g = 1:numel(groups)
    group = groups(g);
    group_path = fullfile(root_dir, group);
    D = dir(group_path);
    subj_folders = D([D.isdir] & ~ismember({D.name},{'.','..','structural'}));

    for s = 1:numel(subj_folders)
        subj = subj_folders(s).name;
        key  = normalize_subject_name(subj);
        subj_path = fullfile(group_path, subj);

        % Alpha
        alpha_file = fullfile(subj_path,"Alpha_estimate.mat");
        if exist(alpha_file,"file")
            A = load(alpha_file);
            alpha_struct = struct( ...
                "Power",   R * log10(10^(-4)+A.Power(:)), ...
                "Width",   R * A.Width(:), ...
                "Exponent",R * A.Exponent(:), ...
                "PAF",     R * A.PAF(:));
            if strcmpi(group,"session1"), alpha_pre_map(key)=alpha_struct;
            else, alpha_post_map(key)=alpha_struct; end
        end

        % Xi
        xi_file = fullfile(subj_path,"Xi_estimate.mat");
        if exist(xi_file,"file")
            X = load(xi_file);
            xi_struct = struct( ...
                "Power",   R * log10(10^(-4)+X.Power(:)), ...
                "Width",   R * X.Width(:), ...
                "Exponent",R * X.Exponent(:));
            if strcmpi(group,"session1"), xi_pre_map(key)=xi_struct;
            else, xi_post_map(key)=xi_struct; end
        end
    end
end

%% Build aligned ROI matrices
Alpha=struct(); Xi=struct();
alpha_fields={"Power","Width","Exponent","PAF"};
xi_fields={"Power","Width","Exponent"};

for k=1:numel(alpha_fields)
    f=alpha_fields{k};
    [Xp,Yp,subs]=build_aligned(alpha_pre_map,alpha_post_map,f);
    if ~isempty(Xp), Alpha.(f).Pre=Xp; Alpha.(f).Post=Yp; Alpha.(f).Subs=subs; end
end
for k=1:numel(xi_fields)
    f=xi_fields{k};
    [Xp,Yp,subs]=build_aligned(xi_pre_map,xi_post_map,f);
    if ~isempty(Xp), Xi.(f).Pre=Xp; Xi.(f).Post=Yp; Xi.(f).Subs=subs; end
end

%% Analysis with ICC(3,1) + bootstrap CI + rigorous permutation test + FDR
icc_results=struct();
processes={"Alpha","Xi"};

for P=1:numel(processes)
    proc=processes{P};
    if strcmp(proc,"Alpha"), params=fieldnames(Alpha); else, params=fieldnames(Xi); end

    for pp=1:numel(params)
        param=params{pp};
        if strcmp(proc,"Alpha")
            Xpre=Alpha.(param).Pre';   % subjects × ROIs
            Xpost=Alpha.(param).Post'; % subjects × ROIs
        else
            Xpre=Xi.(param).Pre';
            Xpost=Xi.(param).Post';
        end

        nsub = size(Xpre,1);
        nrois= size(Xpre,2);

        icc_vals = nan(1,nrois);
        ci_lower = nan(1,nrois);
        ci_upper = nan(1,nrois);
        p_vals   = nan(1,nrois);

        for r = 1:nrois
            data = [Xpre(:,r), Xpost(:,r)];

            % --- Compute observed ICC and bootstrap CI
            [icc_obs, ci] = compute_icc3_1_boot(data, B_boot);
            icc_vals(r) = icc_obs;
            ci_lower(r) = ci(1);
            ci_upper(r) = ci(2);

            % --- Permutation null
            perm_icc = nan(B_perm,1);
            for b=1:B_perm
                perm_idx = randperm(nsub);
                perm_data = [data(:,1), data(perm_idx,2)];
                perm_icc(b) = compute_icc3_1(perm_data);
            end

            % --- Rigorous test: shift null to H0: ICC <= rho0
            shifted_null = perm_icc + rho0;
            if isnan(icc_obs)
                p_vals(r) = NaN;
            elseif icc_obs <= rho0
                p_vals(r) = 1;
            else
                p_vals(r) = (1 + sum(shifted_null >= icc_obs)) / (B_perm + 1);
            end
        end

        % === Apply Benjamini-Hochberg FDR correction ===
        [~,~,~,qvals] = fdr_bh(p_vals,q_FDR);

        key=sprintf("%s_%s",proc,param);
        icc_results.(key).ROIICC   = icc_vals;
        icc_results.(key).CI_lower = ci_lower;
        icc_results.(key).CI_upper = ci_upper;
        icc_results.(key).pvals    = p_vals;
        icc_results.(key).qvals    = qvals;     % corrected q-values
        icc_results.(key).n        = nsub;
        n_sig_FDR = sum(qvals < q_FDR);

        % Median CI across ROIs
        med_ci_lower = nanmedian(ci_lower);
        med_ci_upper = nanmedian(ci_upper);

        fprintf("%s_%s: median ICC=%.3f, 95%% CI=[%.2f, %.2f], rho0=%.2f, median p=%.3f, FDR q<%.2f: %d/%d ROIs significant\n", ...
            proc, param, ...
            nanmedian(icc_vals), ...
            med_ci_lower, med_ci_upper, ...
            rho0, nanmedian(p_vals), q_FDR, ...
            n_sig_FDR, nrois);
    end
end

disp("Rigorous ROI-level ICC + permutation test + FDR finished.");

%% === Helper: ICC(3,1) with bootstrap CI ===
function [icc_obs, ci] = compute_icc3_1_boot(Y, B)
n = size(Y,1);
icc_obs = compute_icc3_1(Y);
boot_icc = nan(B,1);
for b=1:B
    idx = randi(n, n, 1);
    Yb = Y(idx,:);
    boot_icc(b) = compute_icc3_1(Yb);
end
ci = prctile(boot_icc, [2.5, 97.5]);
end

%% === Helper: compute ICC(3,1) ===
function icc = compute_icc3_1(Y)
[n,k] = size(Y);
if n < 2 || k < 2, icc = NaN; return; end
if all(isnan(Y(:))), icc = NaN; return; end

mean_per_subject = mean(Y,2,'omitnan');
mean_per_session = mean(Y,1,'omitnan');
grand_mean = mean(Y(:),'omitnan');

SS_between  = k * sum((mean_per_subject - grand_mean).^2);
SS_sessions = n * sum((mean_per_session - grand_mean).^2);
SS_total    = sum((Y(:) - grand_mean).^2);
SS_error    = SS_total - SS_between - SS_sessions;

if (n-1)==0 || (n-1)*(k-1)==0, icc = NaN; return; end

MS_between = SS_between / (n-1);
MS_error   = SS_error   / ((n-1)*(k-1));

denom = MS_between + (k-1)*MS_error;
if denom <= 0
    icc = NaN;
else
    icc = (MS_between - MS_error) / denom;
end
end

%% === Helper: FDR Benjamini-Hochberg ===
function [h, crit_p, adj_p, qvals] = fdr_bh(pvals,q)
p = pvals(:);
n = sum(~isnan(p));
[p_sorted, sort_ids] = sort(p);
[~, unsort_ids] = sort(sort_ids);
thresh = (1:n)'/n * q;
below = p_sorted <= thresh;
if any(below)
    crit_p = max(p_sorted(below));
else
    crit_p = 0;
end
h = p_sorted <= crit_p;
adj_p = min(1, cummin((n./(1:n)') .* p_sorted,'reverse'));
qvals = nan(size(p));
qvals(sort_ids) = adj_p;
qvals = reshape(qvals,size(pvals));
end

%% === Helper: normalize subject names ===
function key = normalize_subject_name(name)
key = lower(strtrim(name));
end

%% === Helper: align Pre/Post using INTERSECTION ===
function [X_pre,X_post,subs] = build_aligned(pre_map,post_map,field)
pre_keys  = keys(pre_map);
post_keys = keys(post_map);
all_keys  = intersect(pre_keys,post_keys);
all_keys  = sort(all_keys);
nsub = numel(all_keys);
if nsub==0, X_pre=[]; X_post=[]; subs={}; return; end
vlen = numel(pre_map(all_keys{1}).(field));

X_pre  = nan(vlen,nsub);
X_post = nan(vlen,nsub);
for i=1:nsub
    if isKey(pre_map,all_keys{i})
        X_pre(:,i)  = pre_map(all_keys{i}).(field)(:);
    end
    if isKey(post_map,all_keys{i})
        X_post(:,i) = post_map(all_keys{i}).(field)(:);
    end
end
subs = all_keys;
end

%%
% ============================================================
% Combined ICC(3,1) Figure  Bar Plot + Network Heatmap
% Short-term TestRetest Reliability (N = 60, Eyes Closed, 90 min apart)
% ============================================================

colors = [0.2 0.6 0.8;   % eLORETA
          0.8 0.4 0.2;   % LCMV
          0.6 0.8 0.2;   % Xi-AlphaNET
          0.5 0.5 0.5];  % gray

color_xa   = colors(3,:);  % Xi-AlphaNET (a)
color_foof = colors(1,:);  % FOOOF (?)

fields = fieldnames(icc_results);
nFields = numel(fields);

% ============================================================
% Extract ICC summary statistics
% ============================================================
medianICC  = nan(nFields,1);
ci_lower   = nan(nFields,1);
ci_upper   = nan(nFields,1);

for i = 1:nFields
    f = fields{i};
    medianICC(i) = nanmedian(icc_results.(f).ROIICC);
    ci_lower(i)  = nanmedian(icc_results.(f).CI_lower);
    ci_upper(i)  = nanmedian(icc_results.(f).CI_upper);
end

% ============================================================
% Parameter and field order  a first, then ?
% ============================================================
param_labels = {'A\alpha','B\alpha','E\alpha','F\alpha', ...
                'A\xi','B\xi','E\xi'};
param_fields = {'Alpha_Power','Alpha_Width','Alpha_Exponent','Alpha_PAF', ...
                'Xi_Power','Xi_Width','Xi_Exponent'};

order = 1:numel(param_labels);
medianICC = medianICC(order);
ci_lower  = ci_lower(order);
ci_upper  = ci_upper(order);

err_lower = medianICC - ci_lower;
err_upper = ci_upper - medianICC;
bar_colors = [repmat(color_xa,4,1); repmat(color_foof,3,1)];

% ============================================================
% Create Figure
% ============================================================
figure('Color','w','Position',[200 200 1150 500]);

% ============================================================
% LEFT PANEL  Horizontal Bar Plot (a??)
% ============================================================
subplot(1,2,1); hold on;

for i = 1:length(param_labels)
    barh(i, medianICC(i), 'FaceColor', bar_colors(i,:), 'EdgeColor','none');
end

errorbar(medianICC, 1:length(param_labels), err_lower, err_upper, ...
    'horizontal','k','LineStyle','none','LineWidth',1.2);

set(gca,'YTick',1:length(param_labels), ...
         'YTickLabel',param_labels, ...
         'TickLabelInterpreter','tex', ...
         'YDir','reverse', ...
         'FontWeight','bold','FontSize',18);

xlabel('Median ICC(3,1) with 95% bootstrap CI (q_{FDR}<0.05)', ...
       'FontWeight','bold','FontSize',18, 'Interpreter','tex');
ylabel('Spectral parameter','FontWeight','bold','FontSize',18);

title('Whole-brain ICC(3,1) across spectral parameters', ...
      'FontWeight','bold','FontSize',18, 'Interpreter','tex');

xlim([0 1]); ylim([0.5 length(param_labels)+0.5]);
grid on; box on;
set(gca,'LineWidth',1.2,'GridLineStyle',':','GridColor',[0.8 0.8 0.8]);

% ============================================================
% RIGHT PANEL  Heatmap (Parula, High Contrast + Asterisks)
% ============================================================
network_acr = {'VIS','SMN','DAN','VAN','LIM','FPN','DMN'};  % posterior?anterior
n_networks  = numel(network_acr);
n_params    = numel(param_fields);

icc_matrix = nan(n_params,n_networks);
sig_mask   = false(n_params,n_networks);

for p = 1:n_params
    key = param_fields{p};
    vals  = icc_results.(key).ROIICC;
    qvals = icc_results.(key).qvals;
    for n = 1:n_networks
        idx = [2*(n-1)+1, 2*(n-1)+2]; % left & right ROI per network
        icc_matrix(p,n) = median(vals(idx),'omitnan');
        sig_mask(p,n)   = any(qvals(idx) < 0.05);
    end
end

% Flip vertically to match bar plot order
row_order = flip(1:7);
icc_matrix = icc_matrix(row_order,:);
sig_mask   = sig_mask(row_order,:);

subplot(1,2,2);
imagesc(icc_matrix,[0.4 0.9]); hold on;

% --- Add asterisks for non-significant cells ---
for i = 1:n_params
    for j = 1:n_networks
        if ~sig_mask(i,j)
            text(j, i, '*', 'Color','k', 'FontSize',18, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontWeight','bold');
        end
    end
end

% --- Add black grid lines ---
for x = 0.5:1:n_networks+0.5
    plot([x x],[0.5 n_params+0.5],'w','LineWidth',0.6);
end
for y = 0.5:1:n_params+0.5
    plot([0.5 n_networks+0.5],[y y],'w','LineWidth',0.6);
end

% --- Enhanced contrast ---
colormap(parula(256));
cmap = colormap;
x = linspace(0,1,size(cmap,1)).^0.7;  % gamma correction
cmap = interp1(linspace(0,1,numel(x)), cmap, x);
colormap(cmap);

cb = colorbar;
ylabel(cb,'Median ICC(3,1)','FontWeight','bold','FontSize',18);
caxis([0.4 0.9]);

set(gca,'XTick',1:n_networks,'XTickLabel',network_acr, ...
        'YTick',1:n_params,'YTickLabel',param_labels(row_order), ...
        'TickLabelInterpreter','tex', ...
        'YDir','normal','FontWeight','bold','FontSize',18);

xlabel('Yeo-7 network (*q_{FDR}>0.05)', ...
       'FontWeight','bold','FontSize',18, 'Interpreter','tex');
title('Network-wise ICC(3,1) reliability patterns', ...
      'FontWeight','bold','FontSize',18, 'Interpreter','tex');

box on; axis tight;
set(gca,'LineWidth',1.2);

% ============================================================
% Adjust subplot positions (balanced layout)
% ============================================================
subplot(1,2,1); pos1 = get(gca,'Position');
subplot(1,2,2); pos2 = get(gca,'Position');

pos1(1) = 0.06; pos1(3) = 0.38;
pos2(1) = 0.50; pos2(3) = 0.43;

set(subplot(1,2,1),'Position',pos1);
set(subplot(1,2,2),'Position',pos2);
%%
disp("-->> Starting process");

% Order: all Alpha parameters first, then Xi parameters
params = { ...
    'Alpha_Power', ...
    'Alpha_PAF', ...
    'Alpha_Width', ...
    'Alpha_Exponent', ...
    'Xi_Power', ...
    'Xi_Width', ...
    'Xi_Exponent'};

% === LEFT HEMISPHERE VIEW ===
disp('--- Left hemisphere view ---');
for i = 1:numel(params)
    field = params{i};
    if ~isfield(icc_results, field)
        warning("Missing field: %s", field);
        continue;
    end

    J = R_inv * icc_results.(field).ROIICC';
    J(isnan(J)) = 0;

    guide.Visualization.esi_plot_single;
    title(sprintf('%s  Left Hemisphere', strrep(field,'_',' ')), 'Interpreter','none');
    colormap(parula);
    view(-180, 0);   % true LEFT hemisphere (rotate to the left side)
end

% === RIGHT HEMISPHERE VIEW ===
disp('--- Right hemisphere view ---');
for i = 1:numel(params)
    field = params{i};
    if ~isfield(icc_results, field)
        warning("Missing field: %s", field);
        continue;
    end

    J = R_inv * icc_results.(field).ROIICC';
    J(isnan(J)) = 0;

    guide.Visualization.esi_plot_single;
    title(sprintf('%s  Right Hemisphere', strrep(field,'_',' ')), 'Interpreter','none');
    colormap(parula);
    view(0, 0);    % true RIGHT hemisphere
end

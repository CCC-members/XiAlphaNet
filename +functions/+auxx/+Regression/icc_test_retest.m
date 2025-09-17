clc; clear; close all;

%% CONFIG
root_dir       = "/Users/ronald/Downloads/xialphanet_Solutions";
groups         = ["Pre","Mid","Post"];   % now 3 sessions
alpha_sig      = 0.05;   % significance level
rho0           = 0;      % null hypothesis ICC <= rho0
B              = 100;    % number of permutations (increase for real test)
rng(0);

disp("-->> Starting process");

%% Load template cortex and ROI projection matrix
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit

Cortex   = load("templates/Cortex.mat");
template = load("templates/axes.mat");
colorMap = load("templates/mycolormap_brain_basic_conn.mat");
currentAxes = template.axes;
hemis = []; 
[Cortex, iHideVert] = split_hemisphere(Cortex, hemis);

% ROI projection matrix (Nv × 360)
[R,R_inv]=functions.auxx.ModelVectorization.roi_average_operators(Cortex,10);

%% Initialize subject maps for each session
alpha_maps = struct(); xi_maps = struct();
for g = 1:numel(groups)
    alpha_maps.(groups(g)) = containers.Map('KeyType','char','ValueType','any');
    xi_maps.(groups(g))    = containers.Map('KeyType','char','ValueType','any');
end

%% Load data
for g = 1:numel(groups)
    group = groups(g);
    group_path = fullfile(root_dir, group);
    D = dir(group_path);
    subj_folders = D([D.isdir] & ~ismember({D.name},{'.','..','Structural'}));

    for s = 1:numel(subj_folders)
        subj = subj_folders(s).name;
        key  = normalize_subject_name(subj);
        subj_path = fullfile(group_path, subj);

        % Alpha
        alpha_file = fullfile(subj_path,"Alpha_estimate.mat");
        if exist(alpha_file,"file")
            A = load(alpha_file);
            alpha_struct = struct( ...
                "Power",   R * A.Power(:)/mean(1), ...
                "Width",   R * A.Width(:), ...
                "Exponent",R * A.Exponent(:), ...
                "PAF",     R * A.PAF(:));
            alpha_maps.(group)(key) = alpha_struct;
        end

        % Xi
        xi_file = fullfile(subj_path,"Xi_estimate.mat");
        if exist(xi_file,"file")
            X = load(xi_file);
            xi_struct = struct( ...
                "Power",   R * X.Power(:)/mean(1), ...
                "Width",   R * X.Width(:), ...
                "Exponent",R * X.Exponent(:));
            xi_maps.(group)(key) = xi_struct;
        end
    end
end

%% Build aligned ROI matrices
Alpha=struct(); Xi=struct();
alpha_fields={"Power","Width","Exponent","PAF"};
xi_fields={"Power","Width","Exponent"};

for k=1:numel(alpha_fields)
    f=alpha_fields{k};
    [Xcells,subs]=build_aligned_multi(alpha_maps,groups,f);
    if ~isempty(Xcells), Alpha.(f).Data=Xcells; Alpha.(f).Subs=subs; end
end
for k=1:numel(xi_fields)
    f=xi_fields{k};
    [Xcells,subs]=build_aligned_multi(xi_maps,groups,f);
    if ~isempty(Xcells), Xi.(f).Data=Xcells; Xi.(f).Subs=subs; end
end

disp("Summary of aligned ROI data:");
fieldsA=fieldnames(Alpha);
for i=1:numel(fieldsA)
    sz=size(Alpha.(fieldsA{i}).Data{1});
    fprintf("Alpha_%s: %d ROIs × %d subjects × %d sessions\n",fieldsA{i},sz(1),sz(2),numel(Alpha.(fieldsA{i}).Data));
end
fieldsX=fieldnames(Xi);
for i=1:numel(fieldsX)
    sz=size(Xi.(fieldsX{i}).Data{1});
    fprintf("Xi_%s: %d ROIs × %d subjects × %d sessions\n",fieldsX{i},sz(1),sz(2),numel(Xi.(fieldsX{i}).Data));
end

%% Analysis with ICC(3,1) + permutation test
icc_results=struct();
processes={"Alpha","Xi"};

for P=1:numel(processes)
    proc=processes{P};
    if strcmp(proc,"Alpha"), params=fieldnames(Alpha); else, params=fieldnames(Xi); end

    for pp=1:numel(params)
        param=params{pp};
        if strcmp(proc,"Alpha"), Xcells=Alpha.(param).Data;
        else, Xcells=Xi.(param).Data; end

        nsub = size(Xcells{1},2);
        nrois= size(Xcells{1},1);

        icc_vals = nan(1,nrois);
        p_vals   = nan(1,nrois);

        fprintf("\n=== Running %s_%s ===\n",proc,param);

        % Build data matrix: subjects × (sessions*d)
        Y = [];
        for g=1:numel(groups)
            Y = [Y, Xcells{g}']; % subjects × ROIs per session
        end

        [icc_obs, p_vals_roi, pval_global, T_null] = compute_icc3_1_perm_multi(Y, rho0, B);

        fprintf("Null mean=%.4f, std=%.4f\n", mean(T_null), std(T_null));

        key=sprintf("%s_%s",proc,param);
        icc_results.(key).ROIICC=icc_obs.*(1-p_vals_roi);
        icc_results.(key).pvals=p_vals_roi;
        icc_results.(key).pval_global=pval_global;
        icc_results.(key).T_null=T_null;
        icc_results.(key).n=nsub;
    end
end

disp("ROI-level ICC analysis finished.");

%% === Visualization of ALL ROI ICC maps with fixed [0,1] scale ===
fields = fieldnames(icc_results);
for f = 1:numel(fields)
    icc_vals = icc_results.(fields{f}).ROIICC;
    icc_vals = max(0, min(1, icc_vals(:))); % Clip to [0,1]
    J = R_inv * icc_vals;
    figure('Name',fields{f},'Color','w');
    guide.Visualization.esi_plot_single(J);
    title(fields{f}, 'Interpreter','none');
    caxis([0 1]); colorbar;
end


%% === Helper: ICC(3,1) permutation test with within-subject label shuffling ===
function [icc_obs, pvals_roi, pval_global, T_null] = compute_icc3_1_perm_multi(Y, rho0, B)
    % Y: n × (sessions*d)
    [n,m] = size(Y);
    d = m/3; % for 3 sessions
    k = 3;

    % observed ICCs
    icc_obs = nan(1,d);
    for r = 1:d
        icc_obs(r) = compute_icc3_1(Y(:,[r, d+r, 2*d+r]));
    end
    T_obs = max(icc_obs - rho0);

    % permutation distribution
    T_null   = nan(B,1);
    perm_icc = nan(B,d);

    for b=1:B
        Y_perm = Y;
        for i=1:n
            perm_idx = randperm(k); % shuffle session labels for subject i
            Y_perm(i,:) = reshape( Y(i, reshape(1:m, d, k)')',1,[] ); % restructure
            % apply permutation
            blocks = reshape(Y_perm(i,:), d,k);
            blocks = blocks(:,perm_idx);
            Y_perm(i,:) = blocks(:)';
        end

        for r=1:d
            perm_icc(b,r) = compute_icc3_1(Y_perm(:,[r, d+r, 2*d+r]));
        end
        T_null(b) = mean(perm_icc(b,:) - rho0);
    end

    % global p-value
    pval_global = mean(T_null >= T_obs);

    % ROI-wise adjusted p-values
    pvals_roi = nan(1,d);
    for r=1:d
        T_obs_r = icc_obs(r) - rho0;
        pvals_roi(r) = mean(max(perm_icc - rho0,[],2) >= T_obs_r);
    end
end

%% === Helper: compute ICC(3,1) ===
function icc = compute_icc3_1(Y)
[n,k] = size(Y);
mean_per_subject = mean(Y,2);
mean_per_session = mean(Y,1);
grand_mean = mean(Y(:));

SS_between  = k * sum((mean_per_subject - grand_mean).^2);
SS_sessions = n * sum((mean_per_session - grand_mean).^2);
SS_total    = sum((Y(:) - grand_mean).^2);
SS_error    = SS_total - SS_between - SS_sessions;

MS_between = SS_between / (n-1);
MS_error   = SS_error   / ((n-1)*(k-1));

icc = (MS_between - MS_error) / (MS_between + (k-1)*MS_error);
end

%% === Helper: normalize subject names ===
function key = normalize_subject_name(name)
    key = lower(name);                     
    key = strrep(key,'_',' ');             
    key = regexprep(key,'[^a-z\s]','');    
    key = strtrim(key);                    
    tokens = strsplit(key);
    if numel(tokens) >= 2
        key = strjoin(tokens(1:2),' ');    
    end
end

%% === Helper: align N sessions using INTERSECTION ===
function [X_cells, subs] = build_aligned_multi(maps_struct, groups, field)
    all_keys = keys(maps_struct.(groups(1)));
    for g = 2:numel(groups)
        all_keys = intersect(all_keys, keys(maps_struct.(groups(g))));
    end
    all_keys = sort(all_keys);
    nsub = numel(all_keys);

    if nsub==0, X_cells={}; subs={}; return; end
    vlen = numel(maps_struct.(groups(1))(all_keys{1}).(field));

    X_cells = cell(1,numel(groups));
    for g=1:numel(groups)
        Xg = nan(vlen,nsub);
        for i=1:nsub
            Xg(:,i) = maps_struct.(groups(g))(all_keys{i}).(field)(:);
        end
        X_cells{g} = Xg;
    end
    subs = all_keys;
end

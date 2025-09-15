clc; clear; close all;

%% CONFIG
root_dir       = "/Users/ronald/Downloads/xialphanet_Solutions";
groups         = ["Pre","Post"];
alpha_sig      = 0.05;   % significance level
rho0           = 0;    % null hypothesis ICC <= rho0
B              = 10000;   % number of permutations
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
[R,R_inv]=functions.auxx.ModelVectorization.roi_average_operators(Cortex,13);
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
            if strcmpi(group,"Pre"), alpha_pre_map(key)=alpha_struct;
            else, alpha_post_map(key)=alpha_struct; end
        end

        % Xi
        xi_file = fullfile(subj_path,"Xi_estimate.mat");
        if exist(xi_file,"file")
            X = load(xi_file);
            xi_struct = struct( ...
                "Power",   R * X.Power(:)/mean(1), ...
                "Width",   R * X.Width(:), ...
                "Exponent",R * X.Exponent(:));
            if strcmpi(group,"Pre"), xi_pre_map(key)=xi_struct;
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

disp("Summary of aligned ROI data:");
fieldsA=fieldnames(Alpha);
for i=1:numel(fieldsA)
    sz=size(Alpha.(fieldsA{i}).Pre);
    fprintf("Alpha_%s: %d ROIs × %d subjects\n",fieldsA{i},sz(1),sz(2));
end
fieldsX=fieldnames(Xi);
for i=1:numel(fieldsX)
    sz=size(Xi.(fieldsX{i}).Pre);
    fprintf("Xi_%s: %d ROIs × %d subjects\n",fieldsX{i},sz(1),sz(2));
end

%% Analysis with ICC(3,1) + permutation test
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
        p_vals   = nan(1,nrois);

        % Compute ICC + permutation p-value per ROI
        for r = 1:nrois
            data = [Xpre(:,r), Xpost(:,r)];
            [icc_obs, p_val] = compute_icc3_1_perm_multi(data, rho0, B);
            if p_val <= alpha_sig
                icc_vals(r) = icc_obs;
            else
                %icc_vals(r) = 0; % mínimo si no significativo
            end
            p_vals(r) = p_val;
        end

        key=sprintf("%s_%s",proc,param);
        icc_results.(key).ROIICC=icc_vals;
        icc_results.(key).pvals=p_vals;
        icc_results.(key).n=nsub;

        fprintf("%s_%s: median ROI ICC=%.3f, sig=%d/%d, n=%d\n", ...
            proc,param,nanmedian(icc_vals),sum(p_vals<=alpha_sig),nrois,nsub);
    end
end

disp("ROI-level ICC analysis finished.");


%% === Helper: ICC(3,1) permutation test ===
% === new helper ===
function [icc_obs, pvals_roi, pval_global] = compute_icc3_1_perm_multi(Y, rho0, B)
    % Y: n × 2d (Pre,Post for all ROIs concatenated)
    [n,~] = size(Y);
    d = size(Y,2)/2;

    % observed ICCs per ROI
    icc_obs = nan(1,d);
    for r = 1:d
        icc_obs(r) = compute_icc3_1(Y(:,[r, d+r]));
    end
    T_obs = max(icc_obs - rho0);

    % permutation distribution
    T_null = nan(B,1);
    perm_icc = nan(B,d); % store per-ROI too
    for b=1:B
        flips = rand(n,1) > 0.5; % subject-wise flip
        Y_perm = Y;
        Y_perm(flips,:) = Y_perm(flips,[d+1:end, 1:d]); % swap Pre/Post
        for r = 1:d
            perm_icc(b,r) = compute_icc3_1(Y_perm(:,[r, d+r]));
        end
        T_null(b) = mean(perm_icc(b,:) - rho0);
    end

    % global p-value
    pval_global = mean(T_null >= T_obs);

    % ROI-wise adjusted p-values (maxT method)
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

%% === Helper: align Pre/Post using INTERSECTION ===
function [X_pre,X_post,subs] = build_aligned(pre_map,post_map,field)
    pre_keys  = keys(pre_map);
    post_keys = keys(post_map);
    all_keys  = intersect(pre_keys,post_keys); % only subjects in both
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


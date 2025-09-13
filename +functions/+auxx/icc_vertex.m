clc; clear; close all;

%% CONFIG
root_dir       = "/Users/ronald/Downloads/xialphanet_Solutions";
groups         = ["Pre","Post"];
alpha_sig      = 0.05;   % nivel de significación
rho0           = 0.6;    % hipótesis nula ICC <= 0.6
B              = 10000;   % número de permutaciones
rng(0);

disp("-->> Starting process");

%% Initialize subject maps
alpha_pre_map = containers.Map('KeyType','char','ValueType','any');
alpha_post_map= containers.Map('KeyType','char','ValueType','any');
xi_pre_map    = containers.Map('KeyType','char','ValueType','any');
xi_post_map   = containers.Map('KeyType','char','ValueType','any');

%% Load data (cortical average)
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
                "Power",   mean(A.Power(:)), ...
                "Width",   mean(A.Width(:)), ...
                "Exponent",mean(A.Exponent(:)), ...
                "PAF",     mean(A.PAF(:)));
            if strcmpi(group,"Pre"), alpha_pre_map(key)=alpha_struct;
            else, alpha_post_map(key)=alpha_struct; end
        end

        % Xi
        xi_file = fullfile(subj_path,"Xi_estimate.mat");
        if exist(xi_file,"file")
            X = load(xi_file);
            xi_struct = struct( ...
                "Power",   mean(X.Power(:)), ...
                "Width",   mean(X.Width(:)), ...
                "Exponent",mean(X.Exponent(:)));
            if strcmpi(group,"Pre"), xi_pre_map(key)=xi_struct;
            else, xi_post_map(key)=xi_struct; end
        end
    end
end

%% Build aligned subject averages
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

disp("Summary of aligned subject averages:");
fieldsA=fieldnames(Alpha);
for i=1:numel(fieldsA)
    fprintf("Alpha_%s: %d subjects\n",fieldsA{i},numel(Alpha.(fieldsA{i}).Pre));
end
fieldsX=fieldnames(Xi);
for i=1:numel(fieldsX)
    fprintf("Xi_%s: %d subjects\n",fieldsX{i},numel(Xi.(fieldsX{i}).Pre));
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
            Xpre=Alpha.(param).Pre(:);   % subjects × 1
            Xpost=Alpha.(param).Post(:); % subjects × 1
        else
            Xpre=Xi.(param).Pre(:);
            Xpost=Xi.(param).Post(:);
        end

        % Data matrix (nsub × 2)
        data = [Xpre, Xpost];
        icc_obs = compute_icc3_1(data);

        % Permutation distribution
        perm_stats = nan(B,1);
        for b = 1:B
            perm_idx = randperm(size(data,1));
            Yperm = [data(:,1), data(perm_idx,2)];
            perm_stats(b) = compute_icc3_1(Yperm);
        end

        % p-valor unilateral (H0: ICC <= rho0)
        if icc_obs <= rho0
            p_val = 1; % automáticamente no significativo
        else
            p_val = (1 + sum(perm_stats >= icc_obs | perm_stats >= rho0)) / (B+1);
        end

        % Guardar resultados
        key=sprintf("%s_%s",proc,param);
        icc_results.(key).ICC  = icc_obs;
        icc_results.(key).pval = p_val;
        icc_results.(key).n    = size(data,1);

        fprintf("%s_%s: ICC=%.3f, p=%.4f, n=%d\n", ...
            proc,param,icc_obs,p_val,size(data,1));
    end
end

disp("Cortical-average ICC analysis with permutation finished.");

%% Helper: ICC(3,1)
function icc = compute_icc3_1(Y)
Y = Y(all(~isnan(Y),2),:);
[n,k] = size(Y);
if n < 2 || k < 2, icc=NaN; return; end

mean_per_subject = mean(Y,2);
mean_per_session = mean(Y,1);
grand_mean = mean(Y(:));

SS_between  = k * sum((mean_per_subject - grand_mean).^2);
SS_sessions = n * sum((mean_per_session - grand_mean).^2);
SS_total    = sum((Y(:) - grand_mean).^2);
SS_error    = SS_total - SS_between - SS_sessions;

MS_between = SS_between / (n-1);
MS_error   = SS_error   / ((n-1)*(k-1));

den = MS_between + (k-1)*MS_error;
if den <= 0
    icc = NaN;
else
    icc = (MS_between - MS_error) / den;
end
end

%% Helper: normalize subject names
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

%% Helper: align Pre/Post subject averages
function [X_pre,X_post,subs] = build_aligned(pre_map,post_map,field)
    pre_keys  = keys(pre_map);
    post_keys = keys(post_map);
    all_keys  = intersect(pre_keys,post_keys); % solo sujetos en ambos
    all_keys  = sort(all_keys);
    nsub = numel(all_keys);

    if nsub==0, X_pre=[]; X_post=[]; subs={}; return; end

    X_pre  = nan(1,nsub);
    X_post = nan(1,nsub);

    for i=1:nsub
        if isKey(pre_map,all_keys{i})
            X_pre(i)  = pre_map(all_keys{i}).(field)(:);
        end
        if isKey(post_map,all_keys{i})
            X_post(i) = post_map(all_keys{i}).(field)(:);
        end
    end
    subs = all_keys;
end

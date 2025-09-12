clc; clear; close all;

%% CONFIG
root_dir       = "/home/ronaldo/Documents/neuroepo/xialphanet_Solutions";
groups         = ["Pre","Post"];
B              = 100;    % permutations
n_boot         = 1000;     % bootstraps
rng(0);

%% Helper: normalize subject folder names
function key = normalize_subject_name(name)
    key = lower(name);                     % lowercase
    key = strrep(key,'_',' ');             % underscores ? spaces
    key = regexprep(key,'[^a-z\s]','');    % remove non-letters
    key = strtrim(key);                    % trim spaces
    tokens = strsplit(key);
    if numel(tokens) >= 2
        key = strjoin(tokens(1:2),' ');    % keep first 2 words
    end
end

%% Initialize maps
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
                "Power", A.Power(:), ...
                "Width", A.Width(:), ...
                "Exponent", A.Exponent(:), ...
                "PAF", A.PAF(:));
            if strcmpi(group,"Pre"), alpha_pre_map(key)=alpha_struct;
            else, alpha_post_map(key)=alpha_struct; end
        end

        % Xi
        xi_file = fullfile(subj_path,"Xi_estimate.mat");
        if exist(xi_file,"file")
            X = load(xi_file);
            xi_struct = struct( ...
                "Power", X.Power(:), ...
                "Width", X.Width(:), ...
                "Exponent", X.Exponent(:));
            if strcmpi(group,"Pre"), xi_pre_map(key)=xi_struct;
            else, xi_post_map(key)=xi_struct; end
        end
    end
end

%% Helper: align Pre/Post using UNION
function [X_pre,X_post,subs] = build_aligned(pre_map,post_map,field)
    pre_keys  = keys(pre_map);
    post_keys = keys(post_map);
    all_keys  = union(pre_keys,post_keys);
    all_keys  = sort(all_keys);
    nsub = numel(all_keys);

    if ~isempty(pre_keys)
        vlen = numel(pre_map(pre_keys{1}).(field));
    else
        vlen = numel(post_map(post_keys{1}).(field));
    end

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

%% Build aligned matrices
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

disp("Summary of aligned data:");
fieldsA=fieldnames(Alpha);
for i=1:numel(fieldsA)
    sz=size(Alpha.(fieldsA{i}).Pre);
    fprintf("Alpha_%s: %d vertices × %d subjects\n",fieldsA{i},sz(1),sz(2));
end
fieldsX=fieldnames(Xi);
for i=1:numel(fieldsX)
    sz=size(Xi.(fieldsX{i}).Pre);
    fprintf("Xi_%s: %d vertices × %d subjects\n",fieldsX{i},sz(1),sz(2));
end

%% RV-coefficient function
function rv = rv_coefficient(X,Y)
    % drop subjects with NaNs in either session
    mask = all(isfinite(X),2) & all(isfinite(Y),2);
    X=X(mask,:); Y=Y(mask,:);
    % center
    X = X - mean(X,1);
    Y = Y - mean(Y,1);
    A = X*X'; B = Y*Y';
    rv = trace(A*B) / sqrt(trace(A*A)*trace(B*B));
end

%% Analysis
rv_results=struct();
processes={"Alpha","Xi"};

for P=1:numel(processes)
    proc=processes{P};
    if strcmp(proc,"Alpha"), params=fieldnames(Alpha); else, params=fieldnames(Xi); end

    for pp=1:numel(params)
        param=params{pp};
        if strcmp(proc,"Alpha")
            Xpre=Alpha.(param).Pre';  % subjects × vertices
            Xpost=Alpha.(param).Post';
        else
            Xpre=Xi.(param).Pre';
            Xpost=Xi.(param).Post';
        end

        % Observed RV
        rv_obs = rv_coefficient(Xpre,Xpost);

        % Permutation test (shuffle subject correspondence)
        perm_stats=zeros(B,1);
        for b=1:B
            perm_idx=randperm(size(Xpost,1));
            Yb=Xpost(perm_idx,:);
            perm_stats(b)=rv_coefficient(Xpre,Yb);
        end
        p_val=(1+sum(perm_stats>=rv_obs))/(1+B);

        % Bootstrap (resample subjects with replacement)
        nsub=size(Xpre,1);
        boot_stats=zeros(n_boot,1);
        for b=1:n_boot
            idx=randsample(nsub,nsub,true);
            boot_stats(b)=rv_coefficient(Xpre(idx,:),Xpost(idx,:));
        end
        CI=prctile(boot_stats,[2.5,97.5]);

        % Save
        key=sprintf("%s_%s",proc,param);
        rv_results.(key).RV=rv_obs;
        rv_results.(key).p_val=p_val;
        rv_results.(key).CI=CI;
        rv_results.(key).n=size(Xpre,1);

        fprintf("%s_%s: RV=%.3f, p=%.4f, CI=[%.3f,%.3f], n=%d\n", ...
            proc,param,rv_obs,p_val,CI(1),CI(2),size(Xpre,1));
    end
end

disp("Analysis finished with RV.");


%%
%% === Visualization of rv_results with group bars ===
% Assume rv_results is already in workspace

% --- Extract data from struct ---
fields = fieldnames(rv_results);
nParams = numel(fields);

labels_raw = cell(1,nParams);
RV      = zeros(1,nParams);
CI_low  = zeros(1,nParams);
CI_high = zeros(1,nParams);
p_vals  = zeros(1,nParams);

for i = 1:nParams
    key = fields{i};
    labels_raw{i}   = key;
    RV(i)       = rv_results.(key).RV;
    CI_low(i)   = rv_results.(key).CI(1);
    CI_high(i)  = rv_results.(key).CI(2);
    p_vals(i)   = rv_results.(key).p_val;
end

% --- Map raw field names to symbolic labels ---
label_map = containers.Map( ...
    {'Alpha_Power','Alpha_Width','Alpha_Exponent','Alpha_PAF', ...
     'Xi_Power','Xi_Width','Xi_Exponent'}, ...
    {'A\alpha','B\alpha','E\alpha','F\alpha', ...
     'A\xi','B\xi','E\xi'} );

labels = cell(1,nParams);
for i = 1:nParams
    if isKey(label_map, labels_raw{i})
        labels{i} = label_map(labels_raw{i});
    else
        labels{i} = labels_raw{i};
    end
end

% --- Define custom colors ---
alphaColor = [0.2, 0.4, 0.8]; % blue for Alpha
xiColor    = [0.6, 0.8, 0.2]; % green for Xi

colors = zeros(nParams,3);
for i = 1:nParams
    if contains(labels_raw{i},'Alpha')
        colors(i,:) = alphaColor;
    else
        colors(i,:) = xiColor;
    end
end

% --- Figure ---
figure('Color','w','Position',[200 200 900 500]); hold on;
x = 1:nParams;

% Error bars with colored markers
for i = 1:nParams
    errorbar(x(i), RV(i), RV(i)-CI_low(i), CI_high(i)-RV(i), ...
        'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
        'MarkerSize', 8, 'CapSize', 8, 'LineWidth', 1.5);
end

% Significance markers
for i = 1:nParams
    if p_vals(i) < 0.05
        text(x(i), RV(i)+0.07, '*', 'HorizontalAlignment','center', ...
             'Color','r','FontSize',14,'FontWeight','bold');
    else
        text(x(i), RV(i)+0.07, 'ns', 'HorizontalAlignment','center', ...
             'Color',[0.5 0.5 0.5],'FontSize',10,'FontAngle','italic');
    end
end

% Formatting
set(gca,'XTick',x,'XTickLabel',labels,'XTickLabelRotation',0, ...
    'FontSize',13,'FontWeight','bold');
ylabel('RV-coefficient','FontSize',14,'FontWeight','bold');
xlabel('Spectral Parameters (*P<0.05)','FontSize',14,'FontWeight','bold');
ylim([0 1.1]);
xlim([0.5 nParams+0.5]);   % padding left/right
yline(0.5,'--','Color',[0.6 0.6 0.6],'LineWidth',1.2);

title('Test-Retest Reliability (RV-coefficient)', ...
      'FontSize',15,'FontWeight','bold');

% --- Add group bars under x-axis ---
y_bar = -0.08; % vertical position relative to axis (negative = below)
yl = ylim; % store y limits
ylim([yl(1)-0.2 yl(2)]); % extend axis to make space for bars

% Alpha group (positions 1–4)
plot([1 4], [y_bar y_bar], 'k-', 'LineWidth',1.5);
text(mean([1 4]), y_bar-0.05, '\alpha-Process', ...
    'HorizontalAlignment','center','FontSize',13,'FontWeight','bold');

% Xi group (positions 5–7)
plot([5 7], [y_bar y_bar], 'k-', 'LineWidth',1.5);
text(mean([5 7]), y_bar-0.05, '\xi-Process', ...
    'HorizontalAlignment','center','FontSize',13,'FontWeight','bold');

box on; grid on; hold off;

clc; clear; close all;

%% === CONFIGURATION ===
dataDir       = '/Users/ronald/Downloads/FSAve_HCP_MMP1_FSAve_Template_19/data/FSAve_Template/@intra';
mapDir        = '/Users/ronald/Downloads/FSAve_HCP_MMP1_FSAve_Template_19/data/neuromaps';
targetSurface = '@default_subject/tess_cortex_pial_low.mat';

mode      = 'group';   % options: 'group' or 'average'
age_group = 5;           % used only when mode = 'group'
nSpins    = 100;         % number of spins for spin test

%% === LOCATE PARAMETER MAPS ===
switch mode
    case 'group'
        pattern = sprintf('*group%d_sources.mat', age_group);
    case 'average'
        pattern = '*average_sources.mat';
    otherwise
        error('Mode must be ''group'' or ''average''.');
end

files = dir(fullfile(dataDir, pattern));
if isempty(files)
    error('No files found for mode=%s (group=%d) in %s', ...
        mode, age_group, dataDir);
end
fprintf('Found %d parameter maps (%s mode, group=%d)\n', ...
    numel(files), mode, age_group);

%% === RUN NEURORECEPTOR ANALYSIS ===
results_all = struct;
for i = 1:numel(files)
    Yfile = fullfile(files(i).folder, files(i).name);
    fprintf('\n>>> Processing file %d/%d: %s\n', i, numel(files), files(i).name);
    try
        res = functions.auxx.Regression.run_neuroreceptors_analysis( ...
            Yfile, mapDir, targetSurface, nSpins, 'Plot', false);
        results_all(i).file    = files(i).name;
        results_all(i).results = res;
    catch ME
        warning('Error processing %s: %s', files(i).name, ME.message);
        results_all(i).file    = files(i).name;
        results_all(i).results = [];
    end
end

%% === SUMMARY TABLE ===
R2_vals    = arrayfun(@(r) r.results.R2_adj, results_all, 'UniformOutput', true);
pSpin_vals = arrayfun(@(r) r.results.p_spin, results_all, 'UniformOutput', true);

summary = table({results_all.file}', num2cell(R2_vals(:)), num2cell(pSpin_vals(:)), ...
    'VariableNames', {'File','R2_adj','p_spin'});
disp(summary);

%% === DOMINANCE TABLE (CONTRIBUTIONS) ===
M = numel(results_all);  
tags = results_all(1).results.tags; % receptor tags (common to all results)

% --- Normalize tag names to match canonical receptor_order ---
for i = 1:numel(tags)
    switch lower(tags{i})
        case 'a4b2'
            tags{i} = 'α4β2';   % Greek letters
        case 'gabaa'
            tags{i} = 'GABAA';
        case 'vacht'
            tags{i} = 'VAChT';
            % Note: NMDA not present in this dataset, so no mapping needed
    end
end

P = numel(tags);

% --- Build raw contribution matrix ---
contrib_matrix = nan(M,P);
for i = 1:M
    if ~isempty(results_all(i).results)
        contrib_matrix(i,:) = results_all(i).results.contrib(:)' * 100; 
    end
end

% --- Collapse multiple tracers into a single receptor entry ---
unique_tags = unique(tags, 'stable');   
contrib_collapsed = nan(M, numel(unique_tags));
for j = 1:numel(unique_tags)
    cols = strcmp(tags, unique_tags{j});
    contrib_collapsed(:,j) = mean(contrib_matrix(:,cols), 2, 'omitnan');
end
tags = unique_tags;
contrib_matrix = contrib_collapsed;

% --- Define receptor order (excluding NMDA, since it is not in dataset) ---
receptor_order = {'5-HT1a','5-HT1b','5-HT2a','5-HT4','5-HT6','5-HTT', ...
                  'α4β2','CB1','D1','D2','DAT','GABAA','H3','M1', ...
                  'mGluR5','MOR','NET','VAChT'};

% --- Reorder contribution matrix according to receptor_order ---
[~, idx_order] = ismember(receptor_order, tags);
contrib_matrix_full = nan(M, numel(receptor_order));
for j = 1:numel(receptor_order)
    if idx_order(j) > 0
        contrib_matrix_full(:,j) = contrib_matrix(:, idx_order(j));
    end
end
tags_ordered = receptor_order;

% --- Extract R² and p-values ---
R2_adj_vals = R2_vals(:);
pvals       = pSpin_vals(:);

%% === FIGURE: Adjusted R² bars + Dominance heatmap ===
figure('Color','w','Position',[100 100 1200 450]);
tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

% --- Panel 1: Adjusted R² bars ---
nexttile(1);
barh(1:M, R2_adj_vals, 0.4, ...
    'FaceColor',[0.9 0.5 0.1], 'EdgeColor','none');
set(gca,'YDir','reverse','YTick',[]);

row_labels_full = { ...
    '\bf{A\xi}', ...
    '\bf{B\xi}', ...
    '\bf{E\xi}', ...
    '\bf{A\alpha}', ...
    '\bf{B\alpha}', ...
    '\bf{E\alpha}', ...
    '\bf{F\alpha}'};

for i = 1:M
    text(-0.02, i, row_labels_full{i}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'Interpreter','tex','FontSize',12);
end

xlabel('Adjusted R^2   *P_{spin}<0.05','Interpreter','tex');
xlim([0 1]);

ax = gca;
ax.YAxis.TickLength = [0 0];
ax.YAxis.Color = 'w';
ax.XColor = 'k';
box off;
title('Model Fit per Spectral Parameter','FontWeight','bold');

% --- Panel 2: Dominance heatmap ---
nexttile(2);
imagesc(contrib_matrix_full, [0 max(contrib_matrix_full(:),[],'omitnan')]);

% Colormap Hansen-style (white → orange → red)
cmap_orange = [1 1 1; 
               1 0.9 0.6; 
               1 0.7 0.3; 
               1 0.4 0.1; 
               0.7 0 0];
colormap(gca, interp1(linspace(0,1,size(cmap_orange,1)), cmap_orange, linspace(0,1,256)));

cb = colorbar;
cb.Label.String = '% contribution';
cb.Label.FontWeight = 'bold';

set(gca, 'XTick', 1:numel(tags_ordered), 'XTickLabel', tags_ordered, ...
    'XTickLabelRotation', 45, 'YTick', 1:M, 'YTickLabel', []);
ax = gca;
ax.TickLabelInterpreter = 'tex';
axis ij;
box on;
title('Dominance Analysis of Receptors','FontWeight','bold');

% --- Add significance markers ---
for i = 1:M
    if pvals(i) < 0.05
        text(0 - 0.5, i, '*', ...
            'FontSize',14, 'FontWeight','bold', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','middle');
    end
end

% --- Global title ---
sgtitle('Spectral Parameters vs Neuroreceptor Mapping','FontWeight','bold');


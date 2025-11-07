%% === Simulate scalp activity from random ROIs across cortex ===
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit

% --- Load template electrode positions (full set) ---
chanlocs = readlocs('E:\downloads\eeglab2025.0.0\plugins\dipfit\standard_BEM\elec\standard_1005.elc');

% --- ROI labels from cortex atlas ---
labels_rois = {Cortex.Atlas(Cortex.iAtlas).Scouts.Label};
Nr          = size(parameters.Compact_Model.K, 2);  % number of ROIs
K           = parameters.Compact_Model.K;           % (Ne × Nr)

% --- Channel names from your data ---
labels_scalp = data.dnames;   % e.g. {'Fp1','Fp2','F3', ...}
Ne = size(K,1);

% Sanity check
if Ne ~= numel(labels_scalp)
    warning('Mismatch: K has %d electrodes, labels_scalp has %d names.', ...
        Ne, numel(labels_scalp));
end

%% === Pick random ROIs from different lobes ===
roi_groups.occipital = {'V1','V2','V3','V4','V7','V8','LO1','LO2'};
roi_groups.parietal  = {'7','IPS','SPL','Precuneus'};
roi_groups.temporal  = {'STG','MT','MST','TG','TE','PH'};
roi_groups.frontal   = {'10','46','9','44','45','6','FEF','IFG'};

nPerGroup = 2;  % number of ROIs per group

selectedROIs = {};
for fn = fieldnames(roi_groups)'
    group = fn{1};
    keywords = roi_groups.(group);
    idx = [];
    for k = 1:numel(keywords)
        idx = [idx, find(contains(labels_rois, keywords{k}, 'IgnoreCase', true))];
    end
    idx = unique(idx);
    if ~isempty(idx)
        pick = idx(randperm(numel(idx), min(nPerGroup,numel(idx))));
        selectedROIs = [selectedROIs, labels_rois(pick)];
    end
end

%% === Loop and plot ===
nROIs = numel(selectedROIs);
nCols = 4;
nRows = ceil(nROIs/nCols);

figure('Color','w');
for i = 1:nROIs
    targetROI = selectedROIs{i};
    roi_idx   = find(strcmpi(labels_rois, targetROI));

    % delta activation
    x = zeros(Nr,1);
    x(roi_idx) = 1;

    % forward projection
    y = K * x;

    % map onto full chanlocs template
    valuesFull = nan(1, length(chanlocs));
    allLabels  = {chanlocs.labels};
    for j = 1:numel(labels_scalp)
        idx = find(strcmpi(allLabels, labels_scalp{j}));
        if ~isempty(idx)
            valuesFull(idx) = y(j);
        end
    end

    % subplot topoplot
    subplot(nRows, nCols, i);
    topoplot(valuesFull, chanlocs, ...
        'maplimits','maxmin', ...
        'electrodes','labelpoint', ...
        'nosedir','+Y');

    % safe title method (no conflict with "title" shadowing)
    ax = gca;
    set(get(ax,'Title'), 'String', sprintf('%s (ROI #%d)', targetROI, roi_idx), ...
        'Interpreter','none','FontSize',8);
end

sgtitle('Forward projections from random ROIs (Occipital → Frontal)');

%% Load Cortex template
Cortex = load('/Users/ronald/Documents/GitHub/XiAlphaNet/+templates/Cortex.mat');
atlas = Cortex.Atlas(13);   % HCP-MMP1
scouts = atlas.Scouts;

%% Load Python myelin data
data = load('/Users/ronald/Downloads/myelin_hcpmmp1.mat');
labels = string(data.labels);   % Python labels (already Cortex-style)
values = data.values;

%% Initialize vector aligned to Cortex Atlas order
nROIs = numel(scouts);
myelin_values = nan(nROIs,1);
unmatched = {};

%% Match directly (no cleaning needed)
for i = 1:nROIs
    scout_name = string(scouts(i).Label);  % already like "L_V1_ROI L"
    match_idx = find(labels == scout_name);
    if ~isempty(match_idx)
        myelin_values(i) = values(match_idx);
    else
        unmatched{end+1} = scout_name; %#ok<SAGROW>
    end
end

%% Report matching
fprintf("Matched %d / %d regions\n", nROIs - numel(unmatched), nROIs);
fprintf("Unmatched regions: %d\n", numel(unmatched));
if ~isempty(unmatched)
    disp("Example unmatched labels:");
    disp(unmatched(1:min(10,numel(unmatched))))
end

%% Save into Cortex
Cortex.Atlas(13).MyelinValues = myelin_values;
save('/Users/ronald/Downloads/Cortex_with_myelin.mat','Cortex');

%%
% Anchors in myelin units (adjust to your data range)
vals = [0, 1.3, 1.4, 1.6, 1.7];

% Colors: violet/gray → magenta → green → yellow → red
colors = [
     50   40   70   % dark violet at 0
    120   40  120   % magenta at 0.6
     60  150   60   % green at 1.0
    220  220   60   % yellow at 1.3
    220   60   40   % red at 1.7
] / 255;

% Interpolate
cmap = interp1(vals, colors, linspace(min(vals), max(vals), 256), 'linear');

colormap(cmap);
colorbar;
caxis([0 1.7]);  % lock range to your myelin values

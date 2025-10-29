function plot_topo(data)
%% === Imports and template loading ===
import templates.*
import guide.functions.split_hemisphere
import guide.functions.tess_smooth
import guide.functions.tess_hemisplit

% Load template electrode positions
chanlocs = readlocs('E:\downloads\eeglab2025.0.0\plugins\dipfit\standard_BEM\elec\standard_1005.elc');

% --- Channel names and spectra ---
labels_scalp = data.dnames;    % {'Fp1','Fp2','F3',...}
Spec         = data.Spec;      % [nchan × nfreq]
freqs        = data.freq;      % frequency vector

% --- Use only scalp channels (exclude reference) ---
Spec   = Spec(1:length(labels_scalp),:);         % first 18 rows = scalp electrodes
labels = labels_scalp;         % 18 labels
nchan  = size(Spec,1);

%% === Band topoplots ===
bands = struct();
bands.delta = [1 4];
bands.theta = [4 8];
bands.alpha = [8 13];
bands.beta  = [13 30];
bandNames = fieldnames(bands);

bandMaps = nan(length(chanlocs), numel(bandNames));
for b = 1:numel(bandNames)
    fRange = bands.(bandNames{b});
    idxFreq = find(freqs >= fRange(1) & freqs <= fRange(2));
    s = mean(Spec(:, idxFreq), 2);

    valuesFull = nan(1, length(chanlocs));
    allLabels  = {chanlocs.labels};
    for i = 1:length(labels)
        idx = find(strcmpi(allLabels, labels{i}));
        if ~isempty(idx)
            valuesFull(idx) = s(i);
        end
    end
    bandMaps(:, b) = valuesFull(:);
end

%% === Max-alpha channel for spectra highlighting ===
alphaIdx = find(freqs >= 8 & freqs <= 13);
alphaPower = mean(Spec(:, alphaIdx), 2);
[~, maxChanIdx] = max(alphaPower);

cmap = lines(nchan);

%% === Figure 1: Band topoplots + spectra ===
figure('Color','w');

% Row 1: δ, θ, α, β
for b = 1:numel(bandNames)
    subplot(3,4,b);
    topoplot(bandMaps(:,b), chanlocs, ...
        'maplimits','maxmin', ...
        'electrodes','on', ...     % show dots only
        'nosedir','+Y');
    colorbar;
    set(get(gca,'Title'),'String',[bandNames{b} ' band']);
    set(gca,'Color','w');
end

% Row 2: spectra (first 47 bins)
subplot(3,4,5:8); hold on;
for ch = 1:nchan
    if ch == maxChanIdx
        plot(freqs(1:47), log10(Spec(ch,1:47)), 'Color', cmap(ch,:), ...
            'LineWidth', 2, 'DisplayName', [labels{ch} ' *']);
    else
        plot(freqs(1:47), log10(Spec(ch,1:47)), 'Color', cmap(ch,:), ...
            'LineWidth', 1, 'DisplayName', labels{ch});
    end
end
xlabel('Frequency (Hz)');
ylabel('Log_{10} Power');
set(get(gca,'Title'),'String','Spectra (first 47 bins)');
legend('show','Location','northeastoutside');
grid on; set(gca,'Color','w');

% Mark alpha-peak point and annotate electrode name
[~, alphaPeakIdxRel] = max(Spec(maxChanIdx, alphaIdx));
alphaPeakIdx = alphaIdx(alphaPeakIdxRel);
freqStar = freqs(alphaPeakIdx);
powStar  = log10(Spec(maxChanIdx, alphaPeakIdx));
plot(freqStar, powStar, 'p', 'MarkerSize', 12, ...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
text(freqStar, powStar, ['  ' labels{maxChanIdx}], ...
    'Color','k','FontWeight','bold','VerticalAlignment','bottom');

% Row 3: zoomed alpha spectra (6–13 Hz)
subplot(3,4,9:12); hold on;
alphaZoomIdx = find(freqs >= 6 & freqs <= 13);
for ch = 1:nchan
    if ch == maxChanIdx
        plot(freqs(alphaZoomIdx), log10(Spec(ch,alphaZoomIdx)), 'Color', cmap(ch,:), ...
            'LineWidth', 2, 'DisplayName', [labels{ch} ' *']);
    else
        plot(freqs(alphaZoomIdx), log10(Spec(ch,alphaZoomIdx)), 'Color', cmap(ch,:), ...
            'LineWidth', 1, 'DisplayName', labels{ch});
    end
end
xlabel('Frequency (Hz)');
ylabel('Log_{10} Power');
set(get(gca,'Title'),'String','Alpha spectra (6–13 Hz)');
legend('show','Location','northeastoutside');
grid on; set(gca,'Color','w');

% Mark alpha-peak point and annotate electrode name
plot(freqStar, powStar, 'p', 'MarkerSize', 12, ...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
text(freqStar, powStar, ['  ' labels{maxChanIdx}], ...
    'Color','k','FontWeight','bold','VerticalAlignment','bottom');

%% === Figure 2: topoplots per frequency (6–13 Hz) ===
freqRange = 6:13;   % integer Hz steps
nFreqs = numel(freqRange);
figure('Color','w');
for f = 1:nFreqs
    [~, idxF] = min(abs(freqs - freqRange(f)));
    s = Spec(:, idxF);

    valuesFull = nan(1, length(chanlocs));
    allLabels  = {chanlocs.labels};
    for i = 1:length(labels)
        idx = find(strcmpi(allLabels, labels{i}));
        if ~isempty(idx)
            valuesFull(idx) = s(i);
        end
    end

    subplot(2, ceil(nFreqs/2), f);
    topoplot(valuesFull, chanlocs, ...
        'maplimits','maxmin', ...
        'electrodes','on', ...     % dots only, no names
        'nosedir','+Y');
    colorbar;
    set(get(gca,'Title'),'String',sprintf('%d Hz', freqRange(f)));
    set(gca,'Color','w');
end

end

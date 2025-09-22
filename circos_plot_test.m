%% === INPUTS ===
% load('Cortex.mat');   % contains Cortex.Atlas
DAI = abs(Axi_mean).*sign(angle(Axi_mean));            % 360x360 DAI connectivity matrix

% === Select Atlases ===
atlasMMP = Cortex.Atlas(13);   % HCP-MMP1 atlas (360 regions)
atlasYEO = Cortex.Atlas(10);   % Yeo 7-networks atlas (14 regions: L/R split)

nRegions  = length(atlasMMP.Scouts);   % 360
nNetworks = length(atlasYEO.Scouts);   % 14

partName = cell(nRegions,1);   % HCPMMP1 region names
Class    = zeros(nRegions,1);  % network index for each region

%% === Map each MMP1 region to a Yeo network ===
for r = 1:nRegions
    vR = atlasMMP.Scouts(r).Vertices;
    partName{r} = atlasMMP.Scouts(r).Label;

    bestOverlap = 0;
    bestClass   = 1;   % fallback

    for c = 1:nNetworks
        vC = atlasYEO.Scouts(c).Vertices;
        overlap = numel(intersect(vR, vC));

        if overlap > bestOverlap
            bestOverlap = overlap;
            bestClass   = c;
        end
    end
    Class(r) = bestClass;
end

%% === Build className from Yeo atlas labels ===
num2acr = {'VIS','SM','DA','VA','LIM','FP','DMN'};  % Yeo7 canonical acronyms
className = cell(nNetworks,1);

for c = 1:nNetworks
    lbl = atlasYEO.Scouts(c).Label;   % e.g. '7Networks_1 L'
    tok = regexp(lbl,'7Networks_(\d+)\s+([LR])','tokens','once');
    if ~isempty(tok)
        k    = str2double(tok{1});   % 1..7
        hemi = tok{2};               % 'L' or 'R'
        base = num2acr{k};
        className{c} = [base '-' hemi];
    else
        className{c} = lbl;          % fallback
    end
end

%% === Desired order for circos ===
% Left hemi: FP→VIS CCW (frontal→occipital)
% Right hemi: VIS→FP CCW (occipital→frontal)
classNameOrdered = { ...
    'FP-L','DMN-L','VA-L','DA-L','SM-L','LIM-L','VIS-L', ...
    'VIS-R','LIM-R','SM-R','DA-R','VA-R','DMN-R','FP-R'};

% Build old→new lookup
old2new = zeros(1,nNetworks);
for c = 1:numel(classNameOrdered)
    oldIdx = find(strcmpi(className, classNameOrdered{c}),1);
    if ~isempty(oldIdx)
        old2new(oldIdx) = c;
    else
        warning('Class "%s" not found in atlas labels.', classNameOrdered{c});
    end
end

% Map each region’s class
ClassOrdered = old2new(Class);
ClassOrdered(ClassOrdered==0) = 1;  % safety

% === Paso 1: filtrar valores cercanos a cero ===
lowThr = max(abs(DAI(:)));    % máximo valor absoluto permitido (escalar)
DAIfilt = (DAI);
DAIfilt(abs(DAI) >= lowThr) = 0;  % solo conservar los valores "moderados"

% === Paso 1.5: Quitar diagonal explícitamente ===
DAIfilt(1:size(DAIfilt,1)+1:end) = 0;

% === Paso 2: aplicar percentil para sacar ruido muy cerca de 0 ===
vals = nonzeros(DAIfilt(:));   % solo los que sobrevivieron
if ~isempty(vals)
    prc = 95.0;                        % percentil alto
    thrVal = prctile(abs(vals), prc);  % umbral dentro del rango filtrado
    
    DAIfilt(abs(DAIfilt) < thrVal) = 0;
    
    % reforzar: diagonal siempre en cero
    DAIfilt(1:size(DAIfilt,1)+1:end) = 0;
end


% === Colormap (diverging, blue–white–orange) ===
% === Custom diverging colormap: Blue (neg) – Gray (0) – Red (pos) ===
n = 256;
mid = round(n/2);

% Negative side: deep navy → medium blue → gray
neg = [linspace(0,0.2,mid)' ...   % R: almost black
       linspace(0,0.4,mid)' ...   % G: low (keeps blue strong)
       linspace(0.8,1,mid)'];     % B: saturated

% Positive side: gray → bright red
pos = [linspace(0.7,1,mid)' ...   % R: gray → red
       linspace(0.7,0.1,mid)' ... % G: fades out
       linspace(0.7,0.1,mid)'];   % B: fades out

cmap = [neg; pos];

% === Circos plot ===
CC = circosChart(DAIfilt, ClassOrdered, ...
                 'ClassName', classNameOrdered, ...
                 'EdgeColorMap', cmap);
CC = CC.draw();
CC.setClassLabel('FontSize',18,'FontWeight','bold');

% === Node colors (Yeo canonical scheme) ===
colors = [0.3 0.8 0.8;   % FP
          0.98 0.7 0.4;  % DMN
          0.7 0.5 0.9;   % VA
          0.3 0.75 0.4;  % DA
          0.95 0.45 0.45;% SM
          0.95 0.75 0.9; % LIM
          0.58 0.8 0.98];% VIS

networks = {'FP','DMN','VA','DA','SM','LIM','VIS'};
for k = 1:7
    CC.setColor(find(strcmp(classNameOrdered,[networks{k} '-L'])), colors(k,:));
    CC.setColor(find(strcmp(classNameOrdered,[networks{k} '-R'])), colors(k,:));
end

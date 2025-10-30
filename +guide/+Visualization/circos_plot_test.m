%% === INPUTS ===
Cortex = load("+templates/Cortex_with_myelin.mat");   % contains Cortex.Atlas
DAI = target_DAI;            % 360x360 DAI connectivity matrix
 
% === Select Atlases ===
atlasMMP = Cortex.Cortex.Atlas(13);   % HCP-MMP1 atlas (360 regions)
atlasYEO = Cortex.Cortex.Atlas(10);   % Yeo 7-networks atlas (14 regions: L/R split)
nRegions  = length(atlasMMP.Scouts);   % 360
nNetworks = length(atlasYEO.Scouts);   % 14
 
partName = cell(nRegions,1);   % HCPMMP1 region names
Class    = zeros(nRegions,1);  % network index for each region
 
%% === Map each MMP1 region to a Yeo network ===
for r = 1:nRegions
    vR = atlasMMP.Scouts(r).Vertices;
    partName{r} = atlasMMP.Scouts(r).Label;
    bestOverlap = 0;
    bestClass   = 1;
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
num2acr = {'VIS','SM','DAN','VA','LIM','FP','DMN'};  % Yeo7 canonical acronyms
className = cell(nNetworks,1);
 
for c = 1:nNetworks
    lbl = atlasYEO.Scouts(c).Label;
    tok = regexp(lbl,'7Networks_(\d+)\s+([LR])','tokens','once');
    if ~isempty(tok)
        k    = str2double(tok{1});
        hemi = tok{2};
        base = num2acr{k};
        className{c} = [base '-' hemi];
    else
        className{c} = lbl;
    end
end
 
%% === Desired order for circos (posterior?frontal on left; mirror on right) ===
classNameOrdered = { ...
    'FP-L','DMN-L','VA-L','DAN-L','SM-L','LIM-L','VIS-L', ...  % left hemisphere (posterior?frontal)
    'VIS-R','LIM-R','SM-R','DAN-R','VA-R','DMN-R','FP-R'};     % right hemisphere (frontal?posterior)
 
old2new = zeros(1,nNetworks);
for c = 1:numel(classNameOrdered)
    oldIdx = find(strcmpi(className, classNameOrdered{c}),1);
    if ~isempty(oldIdx)
        old2new(oldIdx) = c;
    else
        warning('Class "%s" not found in atlas labels.', classNameOrdered{c});
    end
end
ClassOrdered = old2new(Class);
ClassOrdered(ClassOrdered==0) = 1;
 
%% === Define base gradient (occipital = yellow, frontal = gray) ===
nYeo = 7;
occipitalColor = [1.0 0.9 0.2];
midColor       = [0.9 0.4 0.1];
frontalColor   = [0.15 0.15 0.15];
ColorYeo = interp1([0 0.5 1], ...
    [occipitalColor; midColor; frontalColor], ...
    linspace(0,1,nYeo));
 
%% === Build mirrored hemispheric color order ===
ColorLeft  = flipud(ColorYeo);  % left: FP?VIS = gray?yellow
ColorRight = ColorYeo;          % right: VIS?FP = yellow?gray
ColorOrder = [ColorLeft; ColorRight];
 
%% === Optional small rotation (45° CCW) ===
ColorOrder = circshift(ColorOrder, -rotSteps, 1);
 
%% === Filter DAI ===
DAIfilt = DAI;
DAIfilt(1:size(DAIfilt,1)+1:end) = 0;   % remove diagonal
DAIfilt(isnan(DAIfilt)) = 0;

%% === Filter DAI ===
DAIfilt = DAI;
DAIfilt(1:size(DAIfilt,1)+1:end) = 0;   % remove diagonal
DAIfilt(isnan(DAIfilt)) = 0;

%% === Split by sign (feedforward vs feedback) ===
DAI_pos = (DAIfilt > 0) .* DAIfilt;     % keep positive entries
DAI_neg = (DAIfilt < 0) .* abs(DAIfilt); % keep negative entries as positive magnitudes

% --- Make both symmetric for circos plotting ---
ConnMask_pos = DAI_pos + DAI_pos.';   % combine upper+lower
ConnMask_neg = DAI_neg + DAI_neg.';

% --- Normalize each for visibility ---
if nnz(ConnMask_pos) > 0
    ConnMask_pos = ConnMask_pos ./ max(ConnMask_pos(:));
end
if nnz(ConnMask_neg) > 0
    ConnMask_neg = ConnMask_neg ./ max(ConnMask_neg(:));
end

fprintf('Positive edges: %d | Negative edges: %d\n', ...
    nnz(ConnMask_pos), nnz(ConnMask_neg));


theta = linspace(0,2*pi,400);
r_inner = 1;

%% === POSITIVE (Feedforward) ===
fprintf('Plotting %d positive DAI connections...\n', nnz(ConnMask_pos));
figure('Color','w','Name','Positive DAI (Feedforward)', ...
       'Position',[100 100 1600 1200]);  % large high-res window

% --- Improve on-screen rendering quality ---
set(gcf, 'Renderer', 'opengl');
set(gcf, 'GraphicsSmoothing', 'on');
set(gca, 'FontSmoothing', 'on', ...
         'FontName','Arial', ...
         'FontSize',22, ...
         'FontWeight','bold');

hold on;
fill(r_inner*cos(theta), r_inner*sin(theta), 'k', 'EdgeColor','none');

% --- Draw positive connections ---
CC1 = guide.Visualization.circosChart(ConnMask_pos, ClassOrdered, ...
    'PartName', partName, ...
    'ClassName', classNameOrdered, ...
    'ColorOrder', ColorOrder);
CC1 = CC1.draw();

% --- Label formatting ---
CC1.setPartLabel('Color',[0.1 0.1 0.1], ...
    'FontName','Arial','FontSize',8,'FontWeight','bold', ...
    'Interpreter','none');  % show underscores literally
CC1.setClassLabel('FontName','Arial','FontSize',22,'FontWeight','bold', ...
    'Interpreter','none');
axis equal off;

% --- Ensure all text shows literal underscores ---
set(findall(gcf,'Type','text'),'Interpreter','none');

%% === NEGATIVE (Feedback) ===
fprintf('Plotting %d negative DAI connections...\n', nnz(ConnMask_neg));
figure('Color','w','Name','Negative DAI (Feedback)', ...
       'Position',[200 150 1600 1200]);

% --- Improve on-screen rendering quality ---
set(gcf, 'Renderer', 'opengl');
set(gcf, 'GraphicsSmoothing', 'on');
set(gca, 'FontSmoothing', 'on', ...
         'FontName','Arial', ...
         'FontSize',22, ...
         'FontWeight','bold');

hold on;
fill(r_inner*cos(theta), r_inner*sin(theta), 'k', 'EdgeColor','none');

% --- Draw negative connections ---
CC2 = guide.Visualization.circosChart(ConnMask_neg, ClassOrdered, ...
    'PartName', partName, ...
    'ClassName', classNameOrdered, ...
    'ColorOrder', ColorOrder);
CC2 = CC2.draw();

% --- Label formatting ---
CC2.setPartLabel('Color',[0.1 0.1 0.1], ...
    'FontName','Arial','FontSize',8,'FontWeight','bold', ...
    'Interpreter','none');
CC2.setClassLabel('FontName','Arial','FontSize',22,'FontWeight','bold', ...
    'Interpreter','none');
axis equal off;

% --- Ensure underscores display literally ---
set(findall(gcf,'Type','text'),'Interpreter','none');

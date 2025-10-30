%==========================================================================
%  Xi-AlphaNET: Full Linearized Inverse Resolution Simulation (PSF + CTF)
%  Comparison: Structural Priors ON vs OFF
%==========================================================================

clear; clc; close all;

%% === Imports ===========================================================
import templates.*
import guide.Visualization.*
import guide.functions.split_hemisphere
import functions.auxx.ModelVectorization.*
import functions.auxx.ZeroInflatedModels.*
import functions.auxx.TOperator.*
import functions.auxx.Refine_Solution.*
import functions.auxx.OptimizedOperations.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*

%% === Load Xi-AlphaNET model ============================================
Cortex = load("templates/Cortex_with_myelin.mat");
Cortex = Cortex.Cortex;

parameters = load("/mnt/Develop/Ronaldo/program_working/xialphanet_newresults22/structural/parameters.mat");

dir_data = '/mnt/Develop/Ronaldo/dev/MultinationalNorms';
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
subject_folder = subject_folders(randi(length(subject_folders))).name;
mat_file_path  = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
data_struct    = load(mat_file_path);
Nw = 47;
freq = data_struct.data_struct.freqrange(1:Nw);
parameters.Data.freq = freq;
parameters.Parallel.T = 1;

%% === RUN 1: Priors ON ==================================================
fprintf('\n=== Running with Priors ON  ===\n');
parameters_on = parameters;
parameters_on.Model.C = 0.04 * parameters.Model.C;
parameters_on.Model.D = 1   * parameters.Model.D;

[T, G] = Teval(parameters_on);
[Nw, Nv, Ne] = deal(size(T,3), size(T,2), size(T,1));
rng('default'); n_select = 40;
idx_freq = sort(randsample(Nw, n_select));
R_avg_on = zeros(Nv, Nv);
for k = 1:n_select
    w = idx_freq(k);
    R_avg_on = R_avg_on + G(:,:,w) * T(:,:,w);
end
R_avg_on = abs(R_avg_on / n_select);

% --- Compute cortical distances
[Cortex_, iHideVert] = split_hemisphere(Cortex, []);
Vertices = Cortex_.Vertices; Vertices(iHideVert,:) = [];
distMat = pdist2(Vertices, Vertices);
r_local = 0.07;

% --- Compute PSF metrics
SD_on = zeros(Nv,1); PLE_on = zeros(Nv,1); RSA_on = zeros(Nv,1);
for j = 1:Nv
    psf = R_avg_on(:,j);
    power = abs(psf).^2; power = power / (sum(power)+eps);
    SD_on(j)  = sqrt(sum(distMat(:,j).^2 .* power));
    [~, i_max] = max(abs(psf));
    PLE_on(j) = distMat(i_max, j);
    mask_local = distMat(:,j) <= r_local;
    RSA_on(j)  = 1 - sum(power(mask_local));
end
SD_on = SD_on*1000; PLE_on = PLE_on*1000; RSA_on = RSA_on*100;
fprintf('Priors ON computed successfully.\n');

%% === RUN 2: Priors OFF =================================================
fprintf('\n=== Running with Priors OFF ===\n');
parameters_off = parameters;
parameters_off.Model.C = 1 * parameters.Model.C;
parameters_off.Model.D = 0 * parameters.Model.D;

[T, G] = Teval(parameters_off);
[Nw, Nv, Ne] = deal(size(T,3), size(T,2), size(T,1));
R_avg_off = zeros(Nv, Nv);
for k = 1:n_select
    w = idx_freq(k);
    R_avg_off = R_avg_off + G(:,:,w) * T(:,:,w);
end
R_avg_off = abs(R_avg_off / n_select);

SD_off = zeros(Nv,1); PLE_off = zeros(Nv,1); RSA_off = zeros(Nv,1);
for j = 1:Nv
    psf = R_avg_off(:,j);
    power = abs(psf).^2; power = power / (sum(power)+eps);
    SD_off(j)  = sqrt(sum(distMat(:,j).^2 .* power));
    [~, i_max] = max(abs(psf));
    PLE_off(j) = distMat(i_max, j);
    mask_local = distMat(:,j) <= r_local;
    RSA_off(j)  = 1 - sum(power(mask_local));
end
SD_off = SD_off*1000; PLE_off = PLE_off*1000; RSA_off = RSA_off*100;
fprintf('Priors OFF computed successfully.\n');

%% === Compute summary stats =============================================
data_all = {SD_on, SD_off, PLE_on, PLE_off, RSA_on, RSA_off};
labels_main = {'SD (mm)', 'PLE (mm)', 'RSA (%)'};

for i = 1:6
    [f, x] = ksdensity(data_all{i});
    [~, idx] = max(f);
    modes(i) = x(idx);
    q25(i) = prctile(data_all{i}, 25);
    q75(i) = prctile(data_all{i}, 75);
end

% === Violin Plot (On vs Off grouped) ===================================
colors_on  = [0.64 0.78 0.29];   % green for Prior ON
colors_off = [0.6 0.6 0.6];      % gray for Prior OFF
colors = [colors_on; colors_off; colors_on; colors_off; colors_on; colors_off];

figure('Color','w','Position',[420 240 700 460]);
daviolinplot(data_all, ...
    'violin','half', ...
    'box',3, ...
    'scatter',0, ...
    'violinalpha',0.85, ...
    'colors',colors, ...
    'boxwidth',1.4, 'violinwidth',1.1, ...
    'boxcolors','w', 'boxspacing',1.1);

set(gca,'FontSize',12,'FontWeight','bold','LineWidth',1.4,'TickDir','out');
ylabel('Value','FontSize',14,'FontWeight','bold');
title('Xi-AlphaNET Resolution Metrics: Priors ON vs OFF','FontSize',16,'FontWeight','bold');
grid on; box on;
ylim([0 100]);  % ? Fixed y-axis limit

% --- X-axis labels grouped by metric
x_ticks = [1.5 3.5 5.5];
set(gca,'XTick',x_ticks,'XTickLabel',labels_main);

% --- Top labels "On"/"Off"
label_y = 105;
text(1, label_y, 'On',  'FontSize',13,'FontWeight','bold','Color',colors_on);
text(2, label_y, 'Off', 'FontSize',13,'FontWeight','bold','Color',colors_off);
text(3, label_y, 'On',  'FontSize',13,'FontWeight','bold','Color',colors_on);
text(4, label_y, 'Off', 'FontSize',13,'FontWeight','bold','Color',colors_off);
text(5, label_y, 'On',  'FontSize',13,'FontWeight','bold','Color',colors_on);
text(6, label_y, 'Off', 'FontSize',13,'FontWeight','bold','Color',colors_off);

% --- Add vertical text for mode + interval
y_center = 50;
for i = 1:6
    str = sprintf('Mode = \\bf%.1f\\rm\n[\\bf%.1f-%.1f\\rm]', ...
        modes(i), q25(i), q75(i));
    text(i-0.3, y_center, str, 'Rotation',90, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',11);
end

% --- Legend -------------------------------------------------------------
hold on;
plot(NaN,NaN,'s','MarkerFaceColor',colors_on,'MarkerEdgeColor','none','MarkerSize',10);
plot(NaN,NaN,'s','MarkerFaceColor',colors_off,'MarkerEdgeColor','none','MarkerSize',10);
lgd = legend({'Prior ON','Prior OFF'}, 'FontSize',12, 'Location','northeastoutside');
legend boxoff

fprintf('\nMost probable values (Mode ± IQR):\n');
fprintf('SD  -> ON %.2f [%.2f-%.2f] | OFF %.2f [%.2f-%.2f]\n', ...
    modes(1), q25(1), q75(1), modes(2), q25(2), q75(2));
fprintf('PLE -> ON %.2f [%.2f-%.2f] | OFF %.2f [%.2f-%.2f]\n', ...
    modes(3), q25(3), q75(3), modes(4), q25(4), q75(4));
fprintf('RSA -> ON %.2f [%.2f-%.2f] | OFF %.2f [%.2f-%.2f]\n', ...
    modes(5), q25(5), q75(5), modes(6), q25(6), q75(6));

%%
%% === Cortical Surface Plots: Priors ON vs OFF ==========================
% Common settings
colormap_name = "parula";
scale = [0 100];

% === SD (mm) ============================================================
J =  SD_on; scale = [0 100];
guide.Visualization.esi_plot_hem_R_L;
colormap(colormap_name); caxis(scale);
%colorbar('southoutside','FontSize',11,'FontWeight','bold');
%ylabel(colorbar,'SD (mm) — Prior ON','FontSize',12,'FontWeight','bold');

J = SD_off; scale = [0 100];
guide.Visualization.esi_plot_hem_R_L;
colormap(colormap_name); caxis(scale);
%colorbar('southoutside','FontSize',11,'FontWeight','bold');
%ylabel(colorbar,'SD (mm) — Prior OFF','FontSize',12,'FontWeight','bold');

% === PLE (mm) ===========================================================
J =  PLE_on; scale = [0 100];
guide.Visualization.esi_plot_hem_R_L;
colormap(colormap_name); caxis(scale);
%colorbar('southoutside','FontSize',11,'FontWeight','bold');
%ylabel(colorbar,'PLE (mm) — Prior ON','FontSize',12,'FontWeight','bold');

J =  PLE_off; scale = [0 100];
guide.Visualization.esi_plot_hem_R_L;
colormap(colormap_name); caxis(scale);
%colorbar('southoutside','FontSize',11,'FontWeight','bold');
%ylabel(colorbar,'PLE (mm) — Prior OFF','FontSize',12,'FontWeight','bold');

% === RSA (%) ============================================================
J =  RSA_on; scale = [0 100];
guide.Visualization.esi_plot_hem_R_L;
colormap(colormap_name); caxis(scale);
%colorbar('southoutside','FontSize',11,'FontWeight','bold');
%ylabel(colorbar,'RSA (%) — Prior ON','FontSize',12,'FontWeight','bold');

J =  RSA_off; scale = [0 100];
guide.Visualization.esi_plot_hem_R_L;
colormap(colormap_name); caxis(scale);
%colorbar('southoutside','FontSize',11,'FontWeight','bold');
%ylabel(colorbar,'RSA (%) — Prior OFF','FontSize',12,'FontWeight','bold');

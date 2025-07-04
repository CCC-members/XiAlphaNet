%%
%Svv its the cross
[Nv,~,Nw] = size(Svv);
l = zeros(Nv,Nw);
for j=1:Nw
    l(:,j) = real(diag(Svv(:,:,j)));
end


xi_omega = @(e, omega) e(1) ./ (1 + e(2) * omega.^2).^e(3);
    alpha_omega = @(a, omega) a(1) ./ (1 + a(2) * (omega - a(4)).^2).^a(3);
model_fun_a = @(params, omega) ...
        log(alpha_omega(params(4:7), omega));
model_fun_x = @(params, omega) ...
        log(xi_omega(params(1:3), omega));
model_full = @(params, omega) ...
        log(xi_omega(params(1:3), omega)+alpha_omega(params(4:7), omega));
    % Fit Xi-Alpha model per source
    parfor j = 1:Nv
        j
        %fprintf("Fitting ROI %d/%d...\n", j, Nv);
        try
            parms_j = functions.auxx.GenerateSourceSample.fit_xi_alpha_multi(l(j,:), freq(:), 50, 0);  % doplot = 0
            me(j, :) = parms_j(1:3);   % Xi parameters
            ma(j, :) = parms_j(4:7);   % Alpha parameters
            t_parm{j}  = parms_j;
        catch ME
            warning("Fit failed at ROI %d: %s", j, ME.message);
            me(j, :) = NaN;
            ma(j, :) = NaN;
        end
    end
%%

% Positive colormap: white ? orange ? dark orange
n_gray = 95;     % Gray range from 0 to 0.95
n_color = 5;     % Color range from 0.95 to 1
gray = [0.85, 0.85, 0.85];

clip01 = @(x) min(max(x, 0), 1);

% Colors
dark_orange   = clip01([0.60, 0.30, 0.00]);
bright_orange = clip01([1.00, 0.55, 0.00]);

dark_green    = clip01([0.00, 0.35, 0.10]);
bright_green  = clip01([0.10, 0.85, 0.50]);

dark_blue     = clip01([0.15, 0.30, 0.65]);
bright_blue   = clip01([0.30, 0.60, 1.00]);

% Create gray part (flat gray)
gray_part = repmat(gray, n_gray, 1);

% Positive transitions
orange_positive = [gray_part; ...
                   [linspace(gray(1), bright_orange(1), n_color)', ...
                    linspace(gray(2), bright_orange(2), n_color)', ...
                    linspace(gray(3), bright_orange(3), n_color)']];

green_positive = [gray_part; ...
                  [linspace(gray(1), bright_green(1), n_color)', ...
                   linspace(gray(2), bright_green(2), n_color)', ...
                   linspace(gray(3), bright_green(3), n_color)']];

blue_positive = [gray_part; ...
                 [linspace(gray(1), bright_blue(1), n_color)', ...
                  linspace(gray(2), bright_blue(2), n_color)', ...
                  linspace(gray(3), bright_blue(3), n_color)']];

% Negative transitions (also gray to dark color)
orange_negative = [gray_part; ...
                   [linspace(gray(1), dark_orange(1), n_color)', ...
                    linspace(gray(2), dark_orange(2), n_color)', ...
                    linspace(gray(3), dark_orange(3), n_color)']];

green_negative = [gray_part; ...
                  [linspace(gray(1), dark_green(1), n_color)', ...
                   linspace(gray(2), dark_green(2), n_color)', ...
                   linspace(gray(3), dark_green(3), n_color)']];

blue_negative = [gray_part; ...
                 [linspace(gray(1), dark_blue(1), n_color)', ...
                  linspace(gray(2), dark_blue(2), n_color)', ...
                  linspace(gray(3), dark_blue(3), n_color)']];


% Assumes: freq (1 x Nw), t_parm{1:Nv}, Nv known
Nw = length(freq);
log_alpha_all = nan(Nv, Nw);
log_xi_all    = nan(Nv, Nw);
log_full_all  = nan(Nv, Nw);

% Define model components
xi_omega = @(e, omega) e(1) ./ (1 + e(2) * omega.^2).^e(3);
alpha_omega = @(a, omega) a(1) ./ (1 + a(2) * (omega - a(4)).^2).^a(3);
model_fun_a = @(params, omega) 10*log10(alpha_omega(params(4:7), omega));
model_fun_x = @(params, omega) 10*log10(xi_omega(params(1:3), omega));
model_full  = @(params, omega) 10*log10(xi_omega(params(1:3), omega) + alpha_omega(params(4:7), omega));

% Evaluate models for each source
for j = 1:Nv
    if ~isempty(t_parm{j}) && all(isfinite(t_parm{j}))
        log_alpha_all(j, :) = model_fun_a(t_parm{j}, freq);
        log_xi_all(j, :)    = model_fun_x(t_parm{j}, freq);
        log_full_all(j, :)  = model_full(t_parm{j}, freq);
    end
end

% Compute mean and std across Channels
mean_alpha = nanmean(log_alpha_all, 1);
std_alpha  = nanstd(log_alpha_all, 0, 1);
mean_xi    = nanmean(log_xi_all, 1);
std_xi     = nanstd(log_xi_all, 0, 1);
mean_full  = nanmean(log_full_all, 1);
std_full   = nanstd(log_full_all, 0, 1);

% === Colors ===
orange = orange_negative(end,:);
green  = green_negative(end,:);
black  = [0 0 0];

% === Plot ===
figure('Color','w','Units','normalized','Position',[0.2 0.2 0.6 0.6]);
hold on;

% Shaded error bands
fill([freq, fliplr(freq)], [mean_alpha+std_alpha, fliplr(mean_alpha-std_alpha)], ...
    orange, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([freq, fliplr(freq)], [mean_xi+std_xi, fliplr(mean_xi-std_xi)], ...
    green, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([freq, fliplr(freq)], [mean_full+std_full, fliplr(mean_full-std_full)], ...
    black, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Mean curves
plot(freq, mean_alpha, '-', 'Color', orange, 'LineWidth', 2);
plot(freq, mean_xi,    '-', 'Color', green,  'LineWidth', 2);
plot(freq, mean_full,  '-', 'Color', black,  'LineWidth', 2);
xlim([min(freq) max(freq)])
ylim([5 40])

% Axes and labels
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Power (dB)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

grid on;
box on;

legend({'Alpha ±1σ', 'Xi ±1σ', 'Data ±1σ', ...
        'Alpha Component', 'Xi Component', 'Xi-Alpha Model'}, ...
       'FontSize', 12, 'Location', 'northeast');
title('Average Xi-Alpha Components', 'FontSize', 16, 'FontWeight', 'bold');


% Assumes: freq (1 x Nw), t_parm{1:Nv}, Nv known
Nw = length(freq);
log_alpha_all = nan(Nv, Nw);
log_xi_all    = nan(Nv, Nw);
log_full_all  = nan(Nv, Nw);

% Define model components
xi_omega = @(e, omega) e(1) ./ (1 + e(2) * omega.^2).^e(3);
alpha_omega = @(a, omega) a(1) ./ (1 + a(2) * (omega - a(4)).^2).^a(3);
model_fun_a = @(params, omega) 10*log10(alpha_omega(params(4:7), omega));
model_fun_x = @(params, omega) 10*log10(xi_omega(params(1:3), omega));
model_full  = @(params, omega) 10*log10(xi_omega(params(1:3), omega) + alpha_omega(params(4:7), omega));

% Evaluate models for each source
for j = 1:Nv
    if ~isempty(t_parm{j}) && all(isfinite(t_parm{j}))
        log_alpha_all(j, :) = model_fun_a(t_parm{j}, freq);
        log_xi_all(j, :)    = model_fun_x(t_parm{j}, freq);
        log_full_all(j, :)  = model_full(t_parm{j}, freq);
    end
end

% Compute mean and std across Channels
mean_alpha = nanmean(log_alpha_all, 1);
std_alpha  = nanstd(log_alpha_all, 0, 1);
mean_xi    = nanmean(log_xi_all, 1);
std_xi     = nanstd(log_xi_all, 0, 1);
mean_full  = nanmean(log_full_all, 1);
std_full   = nanstd(log_full_all, 0, 1);

% === Colors ===
orange = orange_negative(end,:);
green  = green_negative(end,:);
black  = [0 0 0];

% === Plot ===
figure('Color','w','Units','normalized','Position',[0.2 0.2 0.6 0.6]);
hold on;

% Shaded error bands
fill([freq, fliplr(freq)], [mean_alpha+std_alpha, fliplr(mean_alpha-std_alpha)], ...
    orange, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([freq, fliplr(freq)], [mean_xi+std_xi, fliplr(mean_xi-std_xi)], ...
    green, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Mean curves
plot(freq, mean_alpha, '-', 'Color', orange, 'LineWidth', 2);
plot(freq, mean_xi,    '-', 'Color', green,  'LineWidth', 2);
xlim([min(freq) max(freq)])
ylim([5 40])
% Axes and labels
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Power (dB)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

grid on;
box on;

legend({'Alpha ±1σ', 'Xi ±1σ', 'Data ±1σ', ...
        'Alpha Component', 'Xi Component', 'Xi-Alpha Model'}, ...
       'FontSize', 12, 'Location', 'northeast');
title('Average Xi-Alpha Components', 'FontSize', 16, 'FontWeight', 'bold');
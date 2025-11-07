function make_receptor_dominance_plot(results_all)
% MAKE_RECEPTOR_DOMINANCE_PLOT
%   Creates summary R² bars and dominance heatmap of receptor contributions.
%
% INPUT
%   results_all : struct array from run_neuroreceptors_analysis
%                 with fields .file and .results

    if isempty(results_all) || isempty(results_all(1).results)
        error('results_all is empty or malformed.');
    end

    % --- Collect values ---
    M      = numel(results_all);   % number of spectral parameters
    tags   = results_all(1).results.tags;
    R2_vals    = arrayfun(@(r) r.results.R2_adj, results_all, 'UniformOutput', true);
    pSpin_vals = arrayfun(@(r) r.results.p_spin, results_all, 'UniformOutput', true);

    % --- Normalize tag names ---
    for i = 1:numel(tags)
        switch lower(tags{i})
            case 'a4b2'
                tags{i} = 'α4β2';
            case 'gabaa'
                tags{i} = 'GABAA';
            case 'vacht'
                tags{i} = 'VAChT';
        end
    end

    % --- Build raw contribution matrix ---
    P = numel(tags);
    contrib_matrix = nan(M,P);
    for i = 1:M
        if ~isempty(results_all(i).results)
            contrib_matrix(i,:) = results_all(i).results.contrib(:)' * 100;
        end
    end

    % --- Collapse duplicate tracers ---
    unique_tags = unique(tags, 'stable');
    contrib_collapsed = nan(M, numel(unique_tags));
    for j = 1:numel(unique_tags)
        cols = strcmp(tags, unique_tags{j});
        contrib_collapsed(:,j) = mean(contrib_matrix(:,cols), 2, 'omitnan');
    end
    tags = unique_tags;
    contrib_matrix = contrib_collapsed;

    % --- Define receptor order (exclude NMDA if not present) ---
    receptor_order = {'5-HT1a','5-HT1b','5-HT2a','5-HT4','5-HT6','5-HTT', ...
                      'α4β2','CB1','D1','D2','DAT','GABAA','H3','M1', ...
                      'mGluR5','MOR','NET','VAChT'};
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

    %% === PLOT ===
    figure('Color','w','Position',[100 100 1200 450]);
    tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

    % Panel 1: Adjusted R² bars
    nexttile(1);
    barh(1:M, R2_adj_vals, 0.6, ...
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
    ax.YAxis.Color = 'w';
    ax.XColor = 'k';
    box off;
    title('Model Fit per Spectral Parameter','FontWeight','bold');

    % Panel 2: Dominance heatmap
    nexttile(2);
    imagesc(contrib_matrix_full, [0 max(contrib_matrix_full(:),[],'omitnan')]);

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
    axis ij; box on;
    title('Dominance Analysis of Receptors','FontWeight','bold');

    % Add significance markers
    for i = 1:M
        if pvals(i) < 0.05
            text(0 - 0.5, i, '*', ...
                'FontSize',14, 'FontWeight','bold', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','middle');
        end
    end

    % Global title
    sgtitle('Spectral Parameters vs Neuroreceptor Mapping','FontWeight','bold');

end

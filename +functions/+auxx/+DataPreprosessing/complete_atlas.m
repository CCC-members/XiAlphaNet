function Atlas = complete_atlas(Cortex, Atlas, mode)
% COMPLETE_ATLAS  Fill uncovered vertices in a Brainstorm atlas.
%
% Atlas = complete_atlas(Cortex, Atlas, mode)
%
% Inputs:
%   Cortex : struct with Cortex.Vertices [Nv x 3] and Cortex.Atlas
%   Atlas  : Atlas struct (e.g. Cortex.Atlas(Cortex.iAtlas))
%   mode   : 'nearest'   = reassign uncovered verts to nearest ROI
%            'unknown'   = add a new ROI called 'Unknown'
%
% Output:
%   Atlas  : modified Atlas with all vertices assigned to an ROI

    Nv = size(Cortex.Vertices,1);

    % === Step 1: mark coverage
    covered = false(Nv,1);
    roi_map = zeros(Nv,1); % ROI index per vertex
    for i = 1:numel(Atlas.Scouts)
        idx = Atlas.Scouts(i).Vertices(:);
        covered(idx) = true;
        roi_map(idx) = i;
    end
    uncovered_idx = find(~covered);

   % fprintf('Atlas covers %d/%d vertices (%.1f%%)\n', ...
    %    sum(covered), Nv, 100*sum(covered)/Nv);
  %  fprintf('Uncovered vertices: %d\n', numel(uncovered_idx));

    if isempty(uncovered_idx)
        disp('No missing vertices. Nothing to do.');
        return;
    end

    % === Step 2: fill uncovered vertices
    switch lower(mode)
        case 'nearest'
            % Build list of covered vertices
            covered_idx = find(covered);

            % Drop any NaN/Inf rows
            validCovered = all(isfinite(Cortex.Vertices(covered_idx,:)),2);
            covered_idx = covered_idx(validCovered);
            validUncovered = all(isfinite(Cortex.Vertices(uncovered_idx,:)),2);
            uncovered_idx = uncovered_idx(validUncovered);

            % Use exhaustive mode (no KDTree, avoids crash)
            idx_nn = knnsearch(Cortex.Vertices(covered_idx,:), ...
                               Cortex.Vertices(uncovered_idx,:), ...
                               'NSMethod','exhaustive');

            for k = 1:numel(uncovered_idx)
                nearest_vert = covered_idx(idx_nn(k));
                nearest_roi  = roi_map(nearest_vert);
                Atlas.Scouts(nearest_roi).Vertices(end+1) = uncovered_idx(k);
                roi_map(uncovered_idx(k)) = nearest_roi;
            end
           % fprintf('All uncovered vertices reassigned to nearest ROI.\n');

        case 'unknown'
            % Create a new ROI for uncovered vertices
            Atlas.Scouts(end+1).Label    = 'Unknown';
            Atlas.Scouts(end).Vertices  = uncovered_idx(:)';
            Atlas.Scouts(end).Seed      = mean(Cortex.Vertices(uncovered_idx,:),1);
            Atlas.Scouts(end).Color     = [0.5 0.5 0.5]; % grey
           % fprintf('Created new ROI "Unknown" with %d vertices.\n', numel(uncovered_idx));

        otherwise
            error('Mode must be ''nearest'' or ''unknown''.');
    end
end

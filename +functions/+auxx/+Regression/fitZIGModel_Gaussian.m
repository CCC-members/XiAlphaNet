function results = fitZIGModel_Gaussian(Y, Age, ROI_labels)
% FITZIGMODEL_GAUSSIAN fits a Zero-Inflated Gaussian model.
%
% Inputs:
%   Y          : [N_voxels x N_subjects] matrix (e.g., PAF or APow)
%   Age        : [N_subjects x 1] vector of subject ages
%   ROI_labels : [N_voxels x 1] numeric array of ROI IDs (optional)
%                If empty, fit single model over all voxels.
%
% Outputs:
%   results : Struct containing ZeroModel and GaussianModel

    if nargin < 3 || isempty(ROI_labels)
        ROI_labels = ones(size(Y,1),1);  % Single group (all voxels together)
    end

    unique_rois = unique(ROI_labels);
    N_rois = length(unique_rois);
    results = struct();

    for r = 1:N_rois
        roi_id = unique_rois(r);
        idx_roi = (ROI_labels == roi_id);
        Y_roi = Y(idx_roi, :);  % Voxels x Subjects

        [V, S] = size(Y_roi);

        % Prepare long format data
        Y_long = Y_roi(:);
        Subject_long = repmat((1:S)', V, 1);
        Age_long = Age(Subject_long)';
        Age2_long = Age_long.^2;

        % ---- 1. Zero-inflation model (logistic regression) ----
        IsZero = (Y_long == 0);

        Zero_tbl = table(Age_long, Age2_long, IsZero);
        ZeroModel = fitglm(Zero_tbl, ...
            'IsZero ~ Age_long + Age2_long', ...
            'Distribution', 'binomial', ...
            'Link', 'logit');

        % ---- 2. Gaussian model for positive values ----
        pos_idx = Y_long > 0;
        Y_pos = Y_long(pos_idx);
        Age_pos = Age_long(pos_idx);
        Age2_pos = Age_pos.^2;

        Pos_tbl = table(Age_pos, Age2_pos, Y_pos);

        GaussianModel = fitglm(Pos_tbl, ...
            'Y_pos ~ Age_pos + Age2_pos', ...
            'Distribution', 'normal');

        % ---- 3. Save Results ----
        roi_fieldname = sprintf('ROI_%d', roi_id);
        results.(roi_fieldname).ZeroModel = ZeroModel;
        results.(roi_fieldname).GaussianModel = GaussianModel;
    end
end

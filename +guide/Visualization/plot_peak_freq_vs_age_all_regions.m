function plot_peak_freq_vs_age_all_regions(dir_results)
% dir_results its the folder that contain the results and the XIALPHANET.json
    % Define the list of regions in the desired order
    regionList = {'Occipital','Parietal', 'Temporal', 'Prefrontal'};
    
    % Define hemispheres with associated colors
    hemispheres = {'L', 'R'};
    hemisphereColors = {'b', 'r'};  % Blue for Left, Red for Right
    
    % Create the figure
    figure('Color', 'w', 'Position', [100, 100, 2000, 1500]);  % Adjust size as needed
    
    % Loop over each region (columns)
    for iReg = 1:length(regionList)
        currentRegion = regionList{iReg};
        
        % Loop over each hemisphere (rows)
        for iHem = 1:length(hemispheres)
            currentHemisphere = hemispheres{iHem};
            currentColor = hemisphereColors{iHem};
            
            % Define row indices
            if iHem == 1  % Left Hemisphere occupies rows 1-3
                rowIndices = 1:3;
            else          % Right Hemisphere occupies rows 4-6
                rowIndices = 4:6;
            end
            
            % Get vertices and main title for the current region and hemisphere
            [vertices, mainTitle] = getVertices(currentRegion, currentHemisphere);
            
            % Collect data for the current region and hemisphere
            [All_Data, all_voxel_activations, all_ages, voxel_positions] = ...
                getDataForRegion(vertices,dir_results);
            
            % Define the age range for evaluation
            age_range = linspace(1, 90, 100)';
            
            % Fit the zero-inflated model
            model_output = fit_zero_inflation_random_effects(all_ages, ...
                                all_voxel_activations, voxel_positions, age_range);
                            
            % Extract predictions
            x         = age_range;
            cond_pred = model_output.predictions.cond;      % E[Y|Y>0]
            cond_se   = model_output.predictions.cond_se;   % Standard Error
            zi_pred   = model_output.predictions.zi;        % P[Y=0]
            zi_se     = model_output.predictions.zi_se;     % Standard Error
            
            % Replace NaNs with mean of non-NaN values
            cond_se(isnan(cond_se)) = nanmean(cond_se(~isnan(cond_se)));
            zi_se(isnan(zi_se))     = nanmean(zi_se(~isnan(zi_se)));
            
            % Handle length discrepancies
            if length(cond_se) ~= length(cond_pred)
                cond_se = min(std(all_voxel_activations - mean(cond_se)), 0.1) * ones(size(x));
                zi_se   = 0.01 * ones(size(x));
            end
            
            % Smooth the standard errors
            cond_se = movmean(cond_se, 20);
            zi_se   = movmean(zi_se, 20);
            
            % Compute final prediction
            final_pred = (1 - zi_pred) .* cond_pred;
            
            % Define scaling factors based on SE mean
            if mean(cond_se) < 1e-2
                sigma1 = 100;
                sigma2 = 10;
            else
                sigma1 = 2;
                sigma2 = 2;
            end
            
            % Calculate upper and lower bounds
            cond_upper = cond_pred + sigma1 * cond_se;
            cond_lower = cond_pred - sigma1 * cond_se;
            
            final_upper = (1 - zi_pred) .* (cond_pred + sigma1 * cond_se);
            final_lower = (1 - zi_pred) .* (cond_pred - sigma1 * cond_se);
            
            zi_final_upper = zi_pred + sigma2 * zi_se;
            zi_final_lower = zi_pred - sigma2 * zi_se;
            
            % Determine subplot position
            % Each region has 6 rows (3 for L, 3 for R)
            % Overall grid is 6 rows x 4 columns
            % For each hemisphere, there are 3 plots stacked vertically
            % Row offset based on hemisphere
            if iHem == 1  % Left Hemisphere
                subplotRowOffset = 0;
            else          % Right Hemisphere
                subplotRowOffset = 3;
            end
            
            % Plot Conditional Prediction (Row 1 for L, Row 4 for R)
            spIndex = subplotRowOffset * length(regionList) + iReg;
            subplot(6, 4, spIndex);
            hold on; grid on;
            fillArea(x, cond_upper, cond_lower, currentColor);
            plot(x, cond_pred, currentColor, 'LineWidth', 2);
            title('E[Y|Y>0]', 'FontWeight', 'bold');
            xlabel('Age (Years)', 'FontWeight', 'bold');
            ylabel('APF (Hz)', 'FontWeight', 'bold');
            ylim([7 9]);
            xlim([0 90])
            %axis tight;
            set(gca, 'FontWeight', 'bold');
            hold off;
            
            % Plot Zero-Inflation Probability (Row 2 for L, Row 5 for R)
            spIndex = (subplotRowOffset + 1) * length(regionList) + iReg;
            subplot(6, 4, spIndex);
            hold on; grid on;
            fillArea(x, zi_final_upper, zi_final_lower, currentColor);
            plot(x, zi_pred, currentColor, 'LineWidth', 2);
            plot(x, 0.5 * ones(size(x)), 'k--', 'LineWidth', 1.5);  % Reference line at 0.5
            title('P[Y=0]', 'FontWeight', 'bold');
            xlabel('Age (Years)', 'FontWeight', 'bold');
            ylabel('Probability', 'FontWeight', 'bold');
            ylim([0 1]);
            xlim([0 90])
            set(gca, 'FontWeight', 'bold');
            hold off;
            
            % Plot Final Prediction (Row 3 for L, Row 6 for R)
            spIndex = (subplotRowOffset + 2) * length(regionList) + iReg;
            subplot(6, 4, spIndex);
            hold on; grid on;
            fillArea(x, final_upper, final_lower, currentColor);
            plot(x, final_pred, currentColor, 'LineWidth', 2);
            title( 'E[Y]', 'FontWeight', 'bold');
            xlabel('Age (Years)', 'FontWeight', 'bold');
            ylabel('APF (Hz)', 'FontWeight', 'bold');
            ylim([2.5 7]);
            xlim([0 90])
            %axis tight;
            set(gca, 'FontWeight', 'bold');
            hold off;
        end
    end
    
    % Adjust subplot spacing
    %tightfig;
    
    % Add a super title for the entire figure
    if exist('sgtitle', 'file')
        sgtitle('Peak Frequency vs Age Across Regions and Hemispheres', 'FontWeight', 'bold', 'FontSize', 16);
    else
        % For older MATLAB versions without sgtitle
        annotation('textbox', [0.5, 0.98, 0, 0], ...
                   'String', 'Peak Frequency vs Age Across Regions and Hemispheres', ...
                   'FontWeight', 'bold', 'FontSize', 16, ...
                   'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
    
%     % Add labels for each column (Regions)
%     for iReg = 1:length(regionList)
%         subplot(6, 4, iReg);
%         if iReg == 1
%             ylabel('Left Hemisphere', 'FontWeight', 'bold', 'Rotation', 90, 'Position', [-10, 0.5, 0]);
%         end
%         subplot(6, 4, (iReg - 1) * 6 + 1);
%         title(regionList{iReg}, 'FontWeight', 'bold', 'FontSize', 14);
%     end
    
    % Add legend or annotations if necessary
end

% --------------------------------------------------------------------
% Helper function to fill the area between upper and lower curves
% --------------------------------------------------------------------
function fillArea(x, upperCurve, lowerCurve, colorName)
    % Ensure column vectors
    x          = x(:);
    upperCurve = upperCurve(:);
    lowerCurve = lowerCurve(:);
    
    % Create filled area
    fill([x; flipud(x)], [upperCurve; flipud(lowerCurve)], colorName, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% --------------------------------------------------------------------
% Helper function to get vertices and main title based on region and hemisphere
% --------------------------------------------------------------------
function [vertices, mainTitle] = getVertices(region, hemisphere)
    % Initialize vertices
    vertices = [];
    
    % Define ROIs based on region and hemisphere
    switch region
        case 'Prefrontal'
            if strcmpi(hemisphere, 'L')
                rois = {'L_10d_ROI L', 'L_10pp_ROI L', 'L_10r_ROI L', 'L_a10p_ROI L'};
            else
                rois = {'R_10d_ROI R', 'R_10pp_ROI R', 'R_10r_ROI R', 'R_a10p_ROI R'};
            end
        case 'Temporal'
            if strcmpi(hemisphere, 'L')
                rois = {'L_TE1a_ROI L', 'L_TE1p_ROI L', 'L_TE2a_ROI L', 'L_TE2p_ROI L'};
            else
                rois = {'R_TE1a_ROI R', 'R_TE1p_ROI R', 'R_TE2a_ROI R', 'R_TE2p_ROI R'};
            end
        case 'Parietal'
            if strcmpi(hemisphere, 'L')
                rois = {'L_7AL_ROI L', 'L_7PC_ROI L', 'L_7PL_ROI L', 'L_7Pm_ROI L'};
            else
                rois = {'R_7AL_ROI R', 'R_7PC_ROI R', 'R_7PL_ROI R', 'R_7Pm_ROI R'};
            end
        case 'Occipital'
            if strcmpi(hemisphere, 'L')
                rois = {'L_V1_ROI L', 'L_V2_ROI L', 'L_V3A_ROI L', 'L_V3B_ROI L'};
            else
                rois = {'R_V1_ROI R', 'R_V2_ROI R', 'R_V3A_ROI R', 'R_V3B_ROI R'};
            end
        otherwise
            error('Invalid region specified. Choose from: Prefrontal, Temporal, Parietal, Occipital.');
    end
    
    % Extract vertices for each ROI
    for i = 1:length(rois)
        vertices = [vertices, getRegionVertices(rois{i})];
    end
    
    % Create a main title string (useful for debugging or logging)
    mainTitle = sprintf('%s (%s)', region, hemisphere);
end

% --------------------------------------------------------------------
% Helper function to collect data for the given set of vertices
% --------------------------------------------------------------------
function [All_Data, all_voxel_activations, all_ages, voxel_positions] = ...
          getDataForRegion(vertices,dir_results)
   
    
    % Load necessary parameters once (if not already loaded)
    
    % Prepare storage
    All_Data = {};
    
    % Initialize index
    index = 1;

   dataset = jsondecode(fileread(strcat(dir_results,'/XIALPHANET.json')));

   for j =1:length(dataset.Participants)
        j;
        participant = dataset.Participants(j);
        participant_age = participant.Age;
        if(isequal(participant.Status,'Completed'))

            %
            All_Data{2,index} =  participant_age;
            Part_Info = jsondecode(fileread(fullfile(dataset.Location,participant.SubID,participant.FileInfo)));
            alpha_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
            a(:,4) = alpha_process.PAF;

            % Apply the calculated threshold to find vertices
            threshold = prctile(a(:,1), 90);
            a(:,4) = a(:,4) .* (a(:,1) > threshold);

            % Extract data from the region of interest
            J = a(vertices, 4);
            J = J(~isnan(J));  % Remove NaNs
        end
        % Store the activation and age
        All_Data{1, index} = J;
        All_Data{2, index} = participant_age;

        index = index + 1;
    end


    % Organize data for the model
    all_voxel_activations = [];
    all_ages             = [];
    voxel_positions      = [];

    for idx = 1:length(All_Data(1, :))
        voxel_data  = All_Data{1, idx};   % Activation for subject
        subject_age = All_Data{2, idx};   % Age
        
        nVox     = numel(voxel_data);
        theseIDs = (1:nVox)';  % Dummy ID for random effect if needed
        
        all_voxel_activations = [all_voxel_activations; voxel_data(:)];
        all_ages              = [all_ages; repmat(subject_age, nVox, 1)];
        voxel_positions       = [voxel_positions; theseIDs];
    end
end

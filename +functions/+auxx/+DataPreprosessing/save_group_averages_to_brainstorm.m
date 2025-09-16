function xialphanet_to_brainstorm(baseDir, inDir, outDir, CortexFile, ageGroups)
% SAVE_XIALPHANET_GROUP_SOURCES
% Convert Xi-AlphaNET group-average parameter maps into Brainstorm-compatible source files.
%
% INPUTS:
%   baseDir     : Base directory for the project (string)
%   inDir       : Directory containing group-average .mat files (string)
%   outDir      : Directory where Brainstorm source files will be saved (string)
%   CortexFile  : Path to cortical surface file (e.g., 'tess_cortex_pial_high_8000V.mat')
%   ageGroups   : Vector of age group indices (e.g., 1:5)
%
% OUTPUTS:
%   Brainstorm source files saved in `outDir` for:
%     - Each parameter × each age group
%     - Each parameter averaged across groups
%
% REQUIREMENTS:
%   - Brainstorm must be on the MATLAB path
%   - Helper function `functions.auxx.DataPreprosessing.save_sources_brainstorm`

% -------------------------------------------------------------------------
% Configuration of processes and parameters
% -------------------------------------------------------------------------
processes.Alpha = {'Power','Bandwidth','Exponent','Frequency'}; % Frequency = PAF
processes.Xi    = {'Power','Bandwidth','Exponent'};             % Xi has no PAF

% Ensure output directory exists
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% -------------------------------------------------------------------------
% Main loop over processes and parameters
% -------------------------------------------------------------------------
procNames = fieldnames(processes);

for ip = 1:numel(procNames)
    proc = procNames{ip};          % 'Alpha' or 'Xi'
    params = processes.(proc);     % List of parameters for this process

    for ipar = 1:numel(params)
        param = params{ipar};

        % Load the .mat file containing all age group variables
        inFile = fullfile(inDir, sprintf('%s_%s.mat', proc, param));
        if ~isfile(inFile)
            warning('Missing file: %s', inFile);
            continue;
        end

        S = load(inFile);  % Load once
        groupMaps = [];    % Container for all groups

        for ig = ageGroups
            varName = sprintf('%s_%s_group%d', proc, param, ig);

            if ~isfield(S, varName)
                warning('Missing variable: %s in %s', varName, inFile);
                continue;
            end

            J = S.(varName);  % Group-level parameter map

            % Save individual group source file
            comment = sprintf('%s %s group%d', proc, param, ig);
            functions.auxx.DataPreprosessing.save_sources_brainstorm( ...
                J, CortexFile, outDir, comment);

            % Collect for group average
            groupMaps(:,ig) = J(:);
        end

        % Save across-group average
        if ~isempty(groupMaps)
            J_avg = mean(groupMaps, 2, 'omitnan');
            comment = sprintf('%s %s average', proc, param);
            functions.auxx.DataPreprosessing.save_sources_brainstorm( ...
                J_avg, CortexFile, outDir, comment);
        end
    end
end

fprintf('[DONE] All source maps saved to: %s\n', outDir);
end

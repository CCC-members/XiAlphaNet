function results_all = run_all_receptor_analysis(dataDir, mapDir, targetSurface, nSpins, mode, age_group, varargin)
% RUN_ALL_RECEPTOR_ANALYSIS
%   Loop through parameter maps (group or average mode), run receptor analysis,
%   and collect results into a struct array.
%
% INPUTS
%   dataDir       : directory containing Xi-AlphaNET parameter .mat maps
%   mapDir        : directory containing receptor maps (neuromaps)
%   targetSurface : Brainstorm surface file (for projection)
%   nSpins        : number of spins for spin test
%   mode          : 'group' or 'average'
%   age_group     : integer, only used if mode='group'
%   varargin      : optional name-value pairs, e.g. 'Plot', false
%
% OUTPUT
%   results_all   : struct array with fields:
%                   .file    (map filename)
%                   .results (struct from run_neuroreceptors_analysis)

    % --- Locate parameter maps ---
    switch mode
        case 'group'
            pattern = sprintf('*group%d_sources.mat', age_group);
        case 'average'
            pattern = '*average_sources.mat';
        otherwise
            error('Mode must be ''group'' or ''average''.');
    end

    files = dir(fullfile(dataDir, pattern));
    if isempty(files)
        error('No files found for mode=%s (group=%d) in %s', ...
            mode, age_group, dataDir);
    end
    fprintf('Found %d parameter maps (%s mode, group=%d)\n', ...
        numel(files), mode, age_group);

    % --- Run analysis ---
    results_all = struct;
    for i = 1:numel(files)
        Yfile = fullfile(files(i).folder, files(i).name);
        fprintf('\n>>> Processing file %d/%d: %s\n', i, numel(files), files(i).name);

        try
            res = run_neuroreceptors_analysis( ...
                Yfile, mapDir, targetSurface, nSpins, varargin{:});
            results_all(i).file    = files(i).name;
            results_all(i).results = res;
        catch ME
            warning('Error processing %s: %s', files(i).name, ME.message);
            results_all(i).file    = files(i).name;
            results_all(i).results = [];
        end
    end
end

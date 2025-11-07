function save_sources_brainstorm(J, CortexFile, outDir, comment)
%SAVE_SOURCES_BRAINSTORM  Save a Brainstorm-compatible sources file
%
%   save_sources_brainstorm(J, CortexFile, outDir, comment)
%
%   INPUTS:
%       J          - Vector (Nv x 1) with values per vertex
%       CortexFile - Full path to the Brainstorm cortex surface (.mat)
%       outDir     - Output directory where results will be saved
%       comment    - Text label shown in Brainstorm (e.g., 'Alpha Power')
%
%   The file is saved in outDir with the name:
%       results_<comment>_sources.mat

    % --- Load cortex to validate dimensions ---
    Cortex = load(CortexFile);
    Nv = size(Cortex.Vertices,1);

    if numel(J) ~= Nv
        error('The vector J must have %d values (one per vertex).', Nv);
    end

    % --- Build Sources structure ---
    Sources = struct();
    Sources.ImagingKernel  = [];
    Sources.ImageGridAmp   = [J(:) J(:)];  % Nv x 2, duplicate columns
    Sources.Std            = [];
    Sources.Whitenor       = [];
    Sources.SourceDecomp   = [];
    Sources.nComponents    = 1;
    Sources.Comment        = comment;
    Sources.Function       = 'sources';
    Sources.DataFile       = [];
    Sources.Time           = [0 1];  % two time samples
    Sources.HeadModelFile  = [];
    Sources.HeadModelType  = 'surface';
    Sources.ChannelFlag    = [];
    Sources.GoodChannel    = [];
    % Use the cortex file path provided
    Sources.SurfaceFile    = CortexFile;
    Sources.Atlas          = [];
    Sources.GridLoc        = Cortex.Vertices;
    Sources.GridOrient     = [];
    Sources.GridAtlas      = [];
    Sources.Options        = [];
    Sources.ColormapType   = [];
    Sources.DisplayUnits   = [];
    Sources.ZScore         = [];
    Sources.nAvg           = 1;
    Sources.History        = {datestr(now), 'save_sources_brainstorm.m', 'Custom import'};
    Sources.Leff           = Nv;

    % --- Ensure output directory exists ---
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % --- Save file in specified output directory ---
    safeComment = regexprep(comment, '\s+', '_'); % replace spaces with underscores
    outFile = fullfile(outDir, ['results_' safeComment '_sources.mat']);
    save(outFile, '-struct', 'Sources');  % save fields directly

    fprintf('Saved Brainstorm sources file: %s\n', outFile);

end

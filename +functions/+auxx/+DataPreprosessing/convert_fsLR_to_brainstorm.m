function convert_fsLR_to_brainstorm()
% CONVERT_FSLR_TO_BRAINSTORM
% Convert fsLR GIFTI (.surf.gii) surfaces into Brainstorm-compatible .mat files.
%
% REQUIREMENTS:
%   - SPM12 installed (provides the "gifti" function)
%   - fsLR surfaces available in neuromaps atlases folder
%
% OUTPUT:
%   Brainstorm-compatible .mat surfaces saved in the specified output folder
%
% Author: Ronald 
% -------------------------------------------------------------------------

    % --- Add SPM12 to path (update if needed) ---
    addpath('/Users/ronald/Downloads/spm');  % path to your SPM12 installation

    % --- Input and output directories ---
    fsLRdir_in  = '/Users/ronald/neuromaps-data/atlases/fsLR';
    fsLRdir_out = '/Users/ronald/Downloads/FSAve_HCP_MMP1_FSAve_Template_19/fsLR';

    if ~exist(fsLRdir_out, 'dir')
        mkdir(fsLRdir_out);
    end

    % --- Which densities to convert ---
    densities = {'4k','32k'};

    % --- Loop over densities and hemispheres ---
    for id = 1:numel(densities)
        den = densities{id};
        for hemi = {'L','R'}
            h = hemi{1};

            % Input fsLR mid-thickness surface
            giiFile = fullfile(fsLRdir_in, ...
                sprintf('tpl-fsLR_den-%s_hemi-%s_midthickness.surf.gii', den, h));

            if ~isfile(giiFile)
                warning('Missing file: %s', giiFile);
                continue;
            end

            % Load surface using gifti (from SPM12)
            g = gifti(giiFile);

            % Build Brainstorm surface structure
            sSurf = struct();
            sSurf.Comment   = sprintf('fsLR-%s %s', den, h);
            sSurf.Vertices  = double(g.vertices);
            sSurf.Faces     = double(g.faces);
            sSurf.FaceVertexCData = [];
            sSurf.Atlas     = [];
            sSurf.Scouts    = [];

            % Output filename
            outFile = fullfile(fsLRdir_out, ...
                sprintf('tess_fsLR_%s_%s.mat', den, h));

            % Save as Brainstorm-compatible .mat
            save(outFile, '-struct', 'sSurf');
            fprintf('[OK] Saved %s\n', outFile);
        end
    end

    fprintf('\n[DONE] fsLR surfaces converted and saved in: %s\n', fsLRdir_out);
end

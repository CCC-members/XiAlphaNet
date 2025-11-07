% ========================================================================
% Full pipeline: Convert fsLR surfaces + import functional map to Brainstorm
% ========================================================================

clc; clear; close all;

%% --- Step 0: Setup ---
% Add SPM to path for gifti support
addpath('/Users/ronald/Downloads/spm');

% Input directories
fsLRdir_in  = '/Users/ronald/neuromaps-data/atlases/fsLR';
mapDir_in   = '/Users/ronald/neuromaps-data/annotations/hcps1200/myelinmap/fsLR';

% Output directory for Brainstorm-ready files
outDir = '/Users/ronald/Downloads/FSAve_HCP_MMP1_FSAve_Template_19/fsLR';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% --- Step 1: Convert hemisphere surfaces (.surf.gii → Brainstorm .mat) ---
dens = '32k';  % choose resolution (4k, 32k, etc.)

for hemi = {'L','R'}
    h = hemi{1};
    giiFile = fullfile(fsLRdir_in, ...
        sprintf('tpl-fsLR_den-%s_hemi-%s_midthickness.surf.gii', dens, h));

    if ~isfile(giiFile)
        error('Missing file: %s', giiFile);
    end

    g = gifti(giiFile);

    sSurf = struct();
    sSurf.Comment  = sprintf('fsLR-%s %s', dens, h);
    sSurf.Vertices = double(g.vertices);
    sSurf.Faces    = double(g.faces);

    outFile = fullfile(outDir, sprintf('tess_fsLR_%s_%s.mat', dens, h));
    save(outFile, '-struct', 'sSurf');
    fprintf('[OK] Saved surface: %s\n', outFile);
end

%% --- Step 2: Merge LH + RH surfaces into full cortex ---
CortexL = load(fullfile(outDir, sprintf('tess_fsLR_%s_L.mat', dens)));
CortexR = load(fullfile(outDir, sprintf('tess_fsLR_%s_R.mat', dens)));

Vertices = [CortexL.Vertices; CortexR.Vertices];
Faces    = [CortexL.Faces; CortexR.Faces + size(CortexL.Vertices,1)];

Cortex = struct();
Cortex.Comment  = sprintf('fsLR-%s Full Cortex', dens);
Cortex.Vertices = Vertices;
Cortex.Faces    = Faces;

CortexFileFull = fullfile(outDir, sprintf('tess_fsLR_%s_full.mat', dens));
save(CortexFileFull, '-struct', 'Cortex');
fprintf('[OK] Merged full cortex saved: %s\n', CortexFileFull);

%% --- Step 3: Load functional map (.func.gii, LH+RH) ---
fileL = fullfile(mapDir_in, ...
    sprintf('source-hcps1200_desc-myelinmap_space-fsLR_den-%s_hemi-L_feature.func.gii', dens));
fileR = fullfile(mapDir_in, ...
    sprintf('source-hcps1200_desc-myelinmap_space-fsLR_den-%s_hemi-R_feature.func.gii', dens));

gL = gifti(fileL);
gR = gifti(fileR);

J = [double(gL.cdata); double(gR.cdata)];
fprintf('[OK] Loaded map data: LH=%d, RH=%d, Total=%d vertices\n', ...
    numel(gL.cdata), numel(gR.cdata), numel(J));

%% --- Step 4: Save map as Brainstorm sources ---
mapOutDir = fullfile(outDir, 'myelinmap');
if ~exist(mapOutDir, 'dir')
    mkdir(mapOutDir);
end

comment = sprintf('HCP-S1200 MyelinMap fsLR-%s', dens);
functions.auxx.DataPreprosessing.save_sources_brainstorm(J, CortexFileFull, mapOutDir, comment);

fprintf('[DONE] Sources saved in: %s\n', mapOutDir);

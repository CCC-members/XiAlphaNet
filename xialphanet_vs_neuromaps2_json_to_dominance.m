clc; clear; close all;

% ========================================================================
%   Xi-AlphaNET vs Neuroreceptors: Full Analysis Pipeline
%
%   End-to-end workflow to link EEG spectral parameters (Xi/AlphaNET)
%   with neuroreceptor PET density maps. The pipeline starts from the
%   raw dataset metadata (.json) and produces final receptor-dominance
%   plots.
%
%   Requirements:
%   --------------------------------------------------------------------
%   (0) A Brainstorm protocol with surface anatomy and cortex meshes
%   (1) Preprocessed EEG data analyzed with Xi-AlphaNET, using the
%       Brainstorm cortex as source space
%   (2) The Neuromaps plugin installed in Brainstorm, with PET
%       neuroreceptor maps downloaded
%
%   Processing steps:
%   --------------------------------------------------------------------
%   (1) Compute group averages of spectral parameters from JSON dataset
%   (2) Convert group averages into Brainstorm-compatible source maps
%   (3) Run receptor regression analysis (Xi/Alpha vs PET densities)
%   (4) Generate dominance plots (Adjusted R² + receptor contributions)
%
%   Author: Ronald García Reyes et al., 2025
% ========================================================================

%% === CONFIGURATION ===
% Base directory for dataset, anatomy, and results
baseDir       = '/Users/ronald/Downloads/FSAve_HCP_MMP1_FSAve_Template_19';

% Input dataset metadata (generated after Xi-AlphaNET processing)
json_path     = fullfile(baseDir, 'XIALPHANET.json');

% Output directories
groupOutDir   = fullfile(baseDir, 'Group_Average_xialphanet'); % Step 1 outputs
bsOutDir      = fullfile(baseDir, 'group_average_sources');    % Step 2 outputs

% Anatomical references (Brainstorm-compatible surfaces)
CortexFile    = fullfile(baseDir, 'anat', 'FSAve_Template', ...
                                      'tess_cortex_pial_high_8000V.mat'); % high-res mesh
targetSurface = '@default_subject/tess_cortex_pial_low.mat';              % low-res mesh for regression

% Neuroreceptor PET maps (provided by Neuromaps plugin)
mapDir        = fullfile(baseDir, 'data', 'neuromaps');

% Spin test configuration
nSpins        = 100;   % number of spin permutations for nonparametric p-values

%% === STEP 1: Compute group averages from JSON ===
% ------------------------------------------------------------------------
% Function: compute_group_averages
% - Loads participant metadata from JSON
% - Extracts Xi and Alpha parameters (Power, Bandwidth, Exponent, Frequency)
% - Stratifies subjects into 5 age groups (0–20, 20–40, 40–60, 60–80, 80–100)
% - Computes group-average maps per spectral parameter
% - Saves results in groupOutDir (e.g., Alpha_Power.mat, Xi_Bandwidth.mat)
% ------------------------------------------------------------------------
fprintf('--- Step 1: Computing group averages from JSON ---\n');
compute_group_averages(json_path, groupOutDir);

%% === STEP 2: Convert group averages into Brainstorm sources ===
% ------------------------------------------------------------------------
% Function: save_group_averages_to_brainstorm
% - Reads group-average parameter maps from Step 1
% - Converts each map into Brainstorm-compatible source structures
% - Saves one .mat per parameter and group
% - Also saves "across-group average" maps for each parameter
% - Output directory: bsOutDir
% ------------------------------------------------------------------------
fprintf('--- Step 2: Saving group averages into Brainstorm format ---\n');
save_group_averages_to_brainstorm(groupOutDir, bsOutDir, CortexFile);

%% === STEP 3: Run receptor regression analysis ===
% ------------------------------------------------------------------------
% Function: run_all_receptor_analysis
% - Loads spectral parameter maps from bsOutDir
% - Loads PET receptor maps from mapDir
% - Performs regression for each parameter against receptor densities
% - Computes Adjusted R², spin-test p-values, and receptor dominance
% - Returns results_all (struct with fits, p-values, receptor weights)
% ------------------------------------------------------------------------
fprintf('--- Step 3: Running receptor regression analysis ---\n');
results_all = run_all_receptor_analysis( ...
    bsOutDir, ...      % directory with Xi/Alpha source maps
    mapDir, ...        % directory with PET receptor maps
    targetSurface, ... % low-res surface for spin tests
    nSpins, ...        % number of spin test permutations
    'average', ...     % mode: 'average' across groups
    1, ...             % group index (not used when mode='average')
    'Plot', false);    % suppress internal plotting

%% === STEP 4: Generate receptor dominance plot ===
% ------------------------------------------------------------------------
% Function: make_receptor_dominance_plot
% - Takes results_all from Step 3
% - Produces a two-panel figure:
%     (1) Bar plot of Adjusted R² per spectral parameter
%     (2) Heatmap of receptor contributions (dominance analysis)
% - Marks parameters with significant fits (* if p_spin < 0.05)
% - Receptor order follows Hansen et al. canonical ordering
% ------------------------------------------------------------------------
fprintf('--- Step 4: Generating dominance plot ---\n');
make_receptor_dominance_plot(results_all);

%% === COMPLETION ===
fprintf('[DONE] Xi-AlphaNET vs Neuroreceptors pipeline completed successfully.\n');

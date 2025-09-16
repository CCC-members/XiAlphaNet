function results = run_neuroreceptors_analysis(Yfile, mapDir, targetSurface, nSpins, varargin)
% RUN_NEURORECEPTORS_ANALYSIS
% Hansen-style receptor mapping with spin test and contribution analysis
%
% INPUTS:
%   Yfile         = Brainstorm .mat file OR numeric vector of vertex values
%   mapDir        = directory with receptor maps (*_sources.mat)
%   targetSurface = Brainstorm surface file (e.g. '@default_subject/tess_cortex_pial_low.mat')
%   nSpins        = number of spins for permutation test
%
% OPTIONAL NAME-VALUE PAIRS:
%   'Method'   : regression method ('ols','robust','ridge') [default: 'ols']
%   'Lambda'   : ridge regression penalty (if Method='ridge') [default: 1]
%   'SaveDir'  : folder to save results & figures [default: [] = no save]
%   'Parallel' : use parallel computing for spins (true/false) [default: false]
%   'Plot'     : show figures (true/false) [default: true]
%
% OUTPUT:
%   results struct with fields:
%       R2_adj, p_spin, beta, contrib, names, tags, R2_perm, config

%% === Parse inputs ===
p = inputParser;
addParameter(p,'Method','ols',@(x)ischar(x)||isstring(x));
addParameter(p,'Lambda',1,@isnumeric);
addParameter(p,'SaveDir',[],@(x)ischar(x)||isstring(x)||isempty(x));
addParameter(p,'Parallel',false,@islogical);
addParameter(p,'Plot',true,@islogical);
parse(p,varargin{:});
method   = lower(p.Results.Method);
lambda   = p.Results.Lambda;
saveDir  = p.Results.SaveDir;
usePar   = p.Results.Parallel;
doPlot   = p.Results.Plot;

%% === Initialize Brainstorm if needed ===
if ~brainstorm('status')
    brainstorm nogui;
end

%% === Load receptor maps ===
maps   = dir(fullfile(mapDir, '**', 'results_surface_*.mat'));
Xfiles = fullfile({maps.folder}, {maps.name});

if isempty(Xfiles)
    error('No receptor maps found in %s', mapDir);
end

fprintf('Found %d receptor maps.\n', numel(Xfiles));

%% === Load or prepare Y data ===
if ischar(Yfile) || isstring(Yfile)
    % Brainstorm results file
    YprojFiles = bst_project_sources({Yfile}, targetSurface, 0, 0);
    sYproj     = in_bst_results(YprojFiles{1}, 1);
    Y          = zscore(sYproj.ImageGridAmp(:,1));
elseif isnumeric(Yfile)
    % Raw vector provided
    Y = zscore(Yfile(:));
else
    error('Yfile must be a Brainstorm file path or numeric vector.');
end
N = length(Y);

%% === Build predictor matrix X ===
P = numel(Xfiles);
X = zeros(N,P);
names = cell(P,1);
tags  = cell(P,1);

for i = 1:P
    tmp = load(Xfiles{i}, 'ImageGridAmp','Comment');
    X(:,i) = zscore(tmp.ImageGridAmp(:,1));

    if isfield(tmp,'Comment')
        names{i} = tmp.Comment;
    else
        [~,names{i}] = fileparts(Xfiles{i});
    end

    % Extract clean tag after "__"
    [~,fname,~] = fileparts(Xfiles{i});
    tokens = regexp(fname, 'results_surface_[^_]+__([^_]+)_', 'tokens');
    if ~isempty(tokens)
        tags{i} = tokens{1}{1};
    else
        tags{i} = fname;
    end
end

%% === Regression ===
switch method
    case 'ols'
        [b,~,~,~,stats] = regress(Y, [ones(N,1), X]);
    case 'robust'
        [b, stats] = robustfit(X, Y);
        b = [b(1); b(2:end)]; % intercept + betas
    case 'ridge'
        b = ridge(Y, X, lambda, 0);
        stats = []; % not directly available
    otherwise
        error('Unknown method: %s', method);
end

% Compute R² and adjusted R²
Yhat = [ones(N,1), X]*b;
SSres = sum((Y-Yhat).^2);
SStot = sum((Y-mean(Y)).^2);
R2    = 1 - SSres/SStot;
R2_adj = 1 - (1-R2)*(N-1)/(N-P-1);
beta   = b(2:end);

%% === Spin test ===
if ischar(Yfile) || isstring(Yfile)
    sSrf   = in_tess_bst(sYproj.SurfaceFile);
    coords = sSrf.Reg.Sphere.Vertices;

    R2_perm = zeros(nSpins,1);

    if usePar
        parfor iSpin = 1:nSpins
            R2_perm(iSpin) = spin_test_single(Y, X, coords, N, P);
        end
    else
        for iSpin = 1:nSpins
            R2_perm(iSpin) = spin_test_single(Y, X, coords, N, P);
        end
    end

    p_spin = (sum(R2_perm >= R2_adj)+1)/(nSpins+1);
else
    warning('Spin test skipped: no surface coords available for raw input.');
    R2_perm = [];
    p_spin = NaN;
end

%% === Contribution analysis ===
if sum(abs(beta)) > 0
    contrib = abs(beta) ./ sum(abs(beta));
else
    contrib = nan(P,1);
end

%% === Pack results ===
results.R2_adj  = R2_adj;
results.p_spin  = p_spin;
results.beta    = beta;
results.contrib = contrib;
results.names   = names;
results.tags    = tags;
results.R2_perm = R2_perm;
results.config  = struct('method',method,'nSpins',nSpins, ...
                         'mapDir',mapDir,'targetSurface',targetSurface);

%% === Display summary ===
fprintf('Adjusted R^2 = %.3f (p_spin = %.3f)\n', R2_adj, p_spin);
disp(table(tags(:), names(:), beta, contrib*100, ...
    'VariableNames', {'Tag','FullName','Beta','PercentContribution'}));

%% === Visualization (optional) ===
if doPlot
    figure;
    bar(R2_adj,'FaceColor',[0.2 0.4 0.8]);
    ylabel('Adjusted R^2'); ylim([0 1]);
    title(sprintf('Model fit (p_{spin}=%.3g)', p_spin));

    figure;
    imagesc(contrib*100);
    colormap(hot); colorbar;
    set(gca,'XTick',1:P,'XTickLabel',tags,'XTickLabelRotation',45);
    ylabel('% Contribution');
    title('Simplified Contribution Analysis');
end

%% === Save results if requested ===
if ~isempty(saveDir)
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    outFile = fullfile(saveDir, sprintf('NeuroreceptorAnalysis_%s.mat', datestr(now,'yyyymmdd_HHMMSS')));
    save(outFile,'results');
    fprintf('Results saved to %s\n', outFile);
end

end

%% === Helper: single spin test iteration ===
function R2p_adj = spin_test_single(Y,X,coords,N,P)
    A = randn(3,3);
    [U,~,V] = svd(A);
    R = U*V';
    if det(R)<0, R(:,1)=-R(:,1); end
    coords_rot = coords*R';
    Mdl = KDTreeSearcher(coords);
    idx = knnsearch(Mdl, coords_rot);
    Y_perm = Y(idx);
    [b,~,~,~,stats_perm] = regress(Y_perm, [ones(N,1), X]);
    R2p = stats_perm(1);
    R2p_adj = 1 - (1-R2p)*(N-1)/(N-P-1);
end

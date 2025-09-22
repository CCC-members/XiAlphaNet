function [data_struct, error_msg] = data_gatherer_flexible(data, SAMPLING_FREQ, cnames, data_code, reference, age, sex, country, eeg_device, keep_signal, properties)
import functions.*
import plugins.*
import plugins.PLGReader.*
import plugins.HarMNqEEG.*

%%% Standard options %%%
desired_order_10_20 = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T3','T4','T5','T6','Fz','Cz','Pz'};
desired_freqres     = 0.390625;   % Hz
desired_epoch_size  = 1./desired_freqres; % seconds
%%%%%%%%%%%%%%%%%%%%%%%%

error_msg   = '';
interpolate = 0;
sp = 1 ./ SAMPLING_FREQ;
[~, epoch_size, nepochs] = size(data);
real_freqres = 1 ./ (sp*epoch_size);

%% --- Adjust epochs if needed for frequency resolution ---
if round(real_freqres*100)/100 < round(desired_freqres*100)/100
    new_epoch_size = floor(desired_epoch_size./sp);
    hm = floor(epoch_size./new_epoch_size);
    ndata = []; p = 1;
    for k=1:nepochs
        for h=1:hm
            ndata(:,:,p) = data(:, (h-1)*new_epoch_size+1:h*new_epoch_size, k);
            p = p+1;
        end
    end
    data = ndata; clear ndata
    [~, epoch_size, nepochs] = size(data);
    real_freqres = 1 ./ (sp*epoch_size);
elseif round(real_freqres*100)/100 > round(desired_freqres*100)/100
    interpolate = 1;
end

%% === Step 1: Channel alignment ===
disp('-->> Aligning channels');

if ~isfield(properties,'channel_params') || ~isfield(properties.channel_params,'labels')
    error_msg = 'properties.channel_params.labels missing';
    data_struct = [];
    return;
end

% Normalize inputs to column cell arrays of char
cnames = cellstr(cnames);  cnames = cnames(:);          % Nx1 cell
labels = cellstr(properties.channel_params.labels);     % target labels
labels = labels(:);                                     % Mx1 cell

% If target is classic 19-ch 10–20, apply alias corrections to INPUT cnames
if numel(labels) == numel(desired_order_10_20)
    disp('-->> Detected 19-channel 10-20 system, applying alias corrections to input names');
    cnames = fix_aliases(cnames);
end

% Align input data to requested labels (case- & space-insensitive)
labels_s = lower(strtrim(string(labels)));
cnames_s = lower(strtrim(string(cnames)));
[found, od] = ismember(labels_s, cnames_s);

if any(~found)
    missing_labels = labels(~found);
    disp('WARNING: Some requested labels not found in input data:');
    disp(missing_labels(:)');
end

% Keep only the labels we can find (order = labels order)
keep_idx = od(found);           % indices into input cnames
data     = data(keep_idx,:,:);  % reorder/keep rows of data
cnames   = labels(found);       % final names (exactly as labels, filtered)

% --- Move reference channel (from properties) to the END (do NOT remove) ---
ref_channel = '';
if isfield(properties.channel_params,'ref_channel')
    ref_channel = properties.channel_params.ref_channel;
end
if ~isempty(ref_channel) && ~strcmpi(ref_channel,'average')
    % Find ref in current cnames (case-insensitive; trim spaces)
    ref_idx = find(strcmpi(strtrim(cnames), strtrim(ref_channel)), 1);
    if ~isempty(ref_idx)
        disp(['-->> Moving reference channel "', cnames{ref_idx}, '" (from properties) to the end']);
        % Build as vertical concatenation (column cells) to avoid horzcat issues
        cnames = [cnames(1:ref_idx-1); cnames(ref_idx+1:end); cnames(ref_idx)];
        data   = cat(1, data([1:ref_idx-1, ref_idx+1:end, ref_idx],:,:));
    else
        disp(['WARNING: Reference channel "', ref_channel, '" (from properties) not found in final channel list']);
    end
end

disp(['-->> Final number of channels: ', num2str(numel(cnames))]);
disp('-->> Final channel order:');
disp(cnames.');

%% === Step 2: Build output structure ===
if nepochs >= 1
    data_struct.name   = data_code;
    data_struct.srate  = SAMPLING_FREQ;
    data_struct.nchan  = numel(cnames);
    data_struct.dnames = cnames;          % column cell
    data_struct.ref    = ref_channel;     % from properties (source of truth)
    data_struct.nt     = epoch_size;
    data_struct.age    = age;
    data_struct.sex    = sex;
    data_struct.pais   = country;
    data_struct.EEGMachine = eeg_device;
    data_struct.nepochs = nepochs;
    if keep_signal
        data_struct.data = data;
    else
        data_struct.data = [];
    end

    %% === Step 3: Compute spectra & cross-spectra ===
    [Spec, data_struct.fmin, data_struct.freqres, data_struct.fmax, ...
        data_struct.CrossM, data_struct.ffteeg] = eeg_to_sp(data, sp);
    data_struct.freqrange = data_struct.fmin:data_struct.freqres:data_struct.fmax;

    %% === Step 4: Interpolation if needed ===
    if interpolate
        % Original (source) frequency grid from eeg_to_sp (DC excluded)
        orig_frange = data_struct.fmin : data_struct.freqres : data_struct.freqres*size(Spec,2);

        % Target (destination) grid: standardized 0.390625 Hz resolution,
        % optionally trimmed to ~19.14 Hz (49*desired_freqres) like your code
        target_frange = desired_freqres : desired_freqres : min(data_struct.fmax, 49*desired_freqres);
        Ntgt = numel(target_frange);

        % --- Interpolate Spec (real, non-negative)
        Spec_i = zeros(size(Spec,1), Ntgt);
        for ch = 1:size(Spec,1)
            Spec_i(ch,:) = pchip(orig_frange, Spec(ch,:), target_frange);
        end

        % Safety check: if interpolation produced negatives, skip subject (your original behavior)
        if any(Spec_i(:) < 0)
            data_struct = [];
            error_msg = 'Interpolation produced negative spectra';
            disp(['SKIPPING subject ', data_code]);
            return;
        end

        % --- Interpolate CrossM (complex): interpolate real and imag parts separately
        nd = size(data_struct.CrossM,1);
        CrossM_i = zeros(nd, nd, Ntgt);
        for i = 1:nd
            for j = 1:nd
                rij = pchip(orig_frange, real(squeeze(data_struct.CrossM(i,j,:))).', target_frange);
                iij = pchip(orig_frange, imag(squeeze(data_struct.CrossM(i,j,:))).', target_frange);
                CrossM_i(i,j,:) = rij + 1i*iij;
            end
        end

        % Commit interpolated data
        data_struct.Spec            = Spec_i;
        data_struct.CrossM          = CrossM_i;
        data_struct.Spec_freqrange  = target_frange;
        data_struct.freqrange       = target_frange;     % <- now identical by construction
        data_struct.freqres         = desired_freqres;
        data_struct.fmin            = target_frange(1);
        data_struct.fmax            = target_frange(end);

    else
        % No interpolation: keep original grid, keep them equal explicitly
        data_struct.Spec            = Spec;
        data_struct.Spec_freqrange  = data_struct.freqrange;  % same grid as CrossM
    end

else
    data_struct = [];
    error_msg = 'No epochs available';
end

end

%% === Helper: alias correction for 10–20 system (robust for cell/string) ===
function cnames = fix_aliases(cnames)
cnames_s = lower(strtrim(string(cnames)));
pos = find(strcmp(cnames_s,'t7'),1); if ~isempty(pos), cnames{pos} = 'T3'; end
pos = find(strcmp(cnames_s,'t8'),1); if ~isempty(pos), cnames{pos} = 'T4'; end
pos = find(strcmp(cnames_s,'p7'),1); if ~isempty(pos), cnames{pos} = 'T5'; end
pos = find(strcmp(cnames_s,'p8'),1); if ~isempty(pos), cnames{pos} = 'T6'; end
end

%% === Spectra / cross-spectra computation (unchanged) ===
function [Spec, fmin, freqres, fmax, CrossM, ffteeg] = eeg_to_sp(EEG, sp)
[nd, nt, ne] = size(EEG);
NFFT   = nt;
freqres = 1 ./ (sp*NFFT);
fmin    = freqres;
fmax    = 1 ./ (2*sp);

ncoefs = floor(nt./2);
ffteeg = zeros(nd, ncoefs, ne);
CrossM = zeros(nd, nd, ncoefs);

for k=1:ne
    A = EEG(:,:,k);
    f = fft(A,NFFT,2);
    f = f(:, 2:ncoefs+1); % drop DC
    ffteeg(:,:,k) = f;
    for h=1:ncoefs
        CrossM(:,:,h) = CrossM(:,:,h) + f(:,h) * f(:,h)';
    end
end

CrossM = CrossM ./ ne;
Spec = zeros(nd, ncoefs);
for k=1:size(CrossM,3)
    Spec(:,k) = diag(CrossM(:,:,k));
end
Spec = real(Spec);
end

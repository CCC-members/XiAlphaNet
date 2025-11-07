function [data_struct, error_msg] = data_gatherer_flexible(data, SAMPLING_FREQ, cnames, data_code, reference, age, sex, country, eeg_device, keep_signal, properties)
import functions.*
import plugins.*
import plugins.PLGReader.*
import plugins.HarMNqEEG.*

%%% Standard options %%%
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

%% --- Step 1: Do NOT subset channels yet. Keep full density ---
disp('-->> Using full channel density for spectra computation');

% Normalize input channel names
cnames = cellstr(cnames);  
cnames = cnames(:);

% Save full channel list before selecting subset later
full_cnames = cnames;

%% --- Step 2: Build output structure ---
if nepochs >= 4
    data_struct.name   = data_code;
    data_struct.srate  = SAMPLING_FREQ;
    data_struct.nchan_full  = numel(full_cnames);  % full density channels
    data_struct.dnames_full = full_cnames;         % store full channel names
    data_struct.ref    = reference;       % just pass through input reference
    data_struct.nt     = epoch_size;
    data_struct.age    = age;
    data_struct.sex    = sex;
    data_struct.pais   = country;
    data_struct.EEGMachine = eeg_device;
    data_struct.nepochs = nepochs;
    if keep_signal
        data_struct.data_full = data;     % keep full signal
    else
        data_struct.data_full = [];
    end

    %% === Step 3: Compute spectra & cross-spectra at full density ===
    [Spec_full, data_struct.fmin, data_struct.freqres, data_struct.fmax, ...
        CrossM_full, ffteeg_full] = eeg_to_sp(data, sp);
    data_struct.freqrange = data_struct.fmin:data_struct.freqres:data_struct.fmax;
    data_struct.Spec_full    = Spec_full;
    data_struct.CrossM_full  = CrossM_full;
    data_struct.ffteeg_full  = ffteeg_full;

    %% === Step 4: Now select target electrodes system ===
    if isfield(properties,'channel_params') && isfield(properties.channel_params,'labels')
        target_labels = cellstr(properties.channel_params.labels);
        target_labels = target_labels(:);

        % case-insensitive matching
        full_lower   = lower(strtrim(string(full_cnames)));
        target_lower = lower(strtrim(string(target_labels)));
        [found, idx] = ismember(target_lower, full_lower);

        if any(~found)
            missing_labels = target_labels(~found);
            disp('WARNING: Some requested labels not found in full channel list:');
            disp(missing_labels(:)');
        end

        idx = idx(found);
        target_labels = target_labels(found);

        % move reference channel last if present
        ref_channel = lower(strtrim(properties.channel_params.ref_channel));
        ref_idx = find(strcmp(lower(strtrim(string(target_labels))), ref_channel), 1);
        if ~isempty(ref_idx)
            order = [setdiff(1:numel(target_labels), ref_idx, 'stable'), ref_idx];
            idx = idx(order);
            target_labels = target_labels(order);
        end

        % subset the full Spec and CrossM
        Spec_sub   = Spec_full(idx,:);
        CrossM_sub = CrossM_full(idx,idx,:);

        data_struct.nchan  = numel(target_labels);
        data_struct.dnames = target_labels;
        data_struct.Spec   = Spec_sub;
        data_struct.CrossM = CrossM_sub;
    else
        % no target channel system specified
        data_struct.nchan  = size(Spec_full,1);
        data_struct.dnames = full_cnames;
        data_struct.Spec   = Spec_full;
        data_struct.CrossM = CrossM_full;
    end

    %% === Step 5: Interpolation if needed (apply to the subset) ===
    if interpolate
        orig_frange = data_struct.fmin : data_struct.freqres : data_struct.freqres*size(Spec_sub,2);
        target_frange = desired_freqres : desired_freqres : min(data_struct.fmax, 49*desired_freqres);
        Ntgt = numel(target_frange);

        Spec_i = zeros(size(data_struct.Spec,1), Ntgt);
        for ch = 1:size(data_struct.Spec,1)
            Spec_i(ch,:) = pchip(orig_frange, data_struct.Spec(ch,:), target_frange);
        end

        nd = size(data_struct.CrossM,1);
        CrossM_i = zeros(nd, nd, Ntgt);
        for i = 1:nd
            for j = 1:nd
                rij = pchip(orig_frange, real(squeeze(data_struct.CrossM(i,j,:))).', target_frange);
                iij = pchip(orig_frange, imag(squeeze(data_struct.CrossM(i,j,:))).', target_frange);
                CrossM_i(i,j,:) = rij + 1i*iij;
            end
        end

        data_struct.Spec            = Spec_i;
        data_struct.CrossM          = CrossM_i;
        data_struct.Spec_freqrange  = target_frange;
        data_struct.freqrange       = target_frange;
        data_struct.freqres         = desired_freqres;
        data_struct.fmin            = target_frange(1);
        data_struct.fmax            = target_frange(end);
    else
        data_struct.Spec_freqrange  = data_struct.freqrange;
    end

else
    data_struct = [];
    error_msg = 'No epochs available';
end

end

%% === Spectra / cross-spectra computation ===
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

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

%% === Step 1: Channel alignment ===
disp('-->> Aligning channels');

if ~isfield(properties,'channel_params') || ~isfield(properties.channel_params,'labels')
    error_msg = 'properties.channel_params.labels missing';
    data_struct = [];
    return;
end

% Normalize inputs to column cell arrays of char
cnames = cellstr(cnames);  cnames = cnames(:);
labels = cellstr(properties.channel_params.labels);
labels = labels(:);

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
keep_idx = od(found);
data     = data(keep_idx,:,:);
cnames   = labels(found);

disp(['-->> Final number of channels: ', num2str(numel(cnames))]);
disp('-->> Final channel order:');
disp(cnames.');

%% === Step 2: Build output structure ===
if nepochs >= 4
    data_struct.name   = data_code;
    data_struct.srate  = SAMPLING_FREQ;
    data_struct.nchan  = numel(cnames);
    data_struct.dnames = cnames;
    data_struct.ref    = reference;       % just pass through input reference
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
        orig_frange = data_struct.fmin : data_struct.freqres : data_struct.freqres*size(Spec,2);
        target_frange = desired_freqres : desired_freqres : min(data_struct.fmax, 49*desired_freqres);
        Ntgt = numel(target_frange);

        Spec_i = zeros(size(Spec,1), Ntgt);
        for ch = 1:size(Spec,1)
            Spec_i(ch,:) = pchip(orig_frange, Spec(ch,:), target_frange);
        end

        if any(Spec_i(:) < 0)
            data_struct = [];
            error_msg = 'Interpolation produced negative spectra';
            disp(['SKIPPING subject ', data_code]);
            return;
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
        data_struct.Spec            = Spec;
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

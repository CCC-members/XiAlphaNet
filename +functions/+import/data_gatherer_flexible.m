function [data_struct, error_msg] = data_gatherer_flexible(data, SAMPLING_FREQ, cnames, data_code, reference, age, sex, country, eeg_device, keep_signal, properties)
import functions.*
import plugins.*
import plugins.PLGReader.*
import plugins.HarMNqEEG.*

%%% Standard options for the MultiNational Norms project %%%
desired_order_10_20 = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T3','T4','T5','T6','Fz','Cz','Pz'};
desired_freqres = 0.390625; % Hz
desired_epoch_size = 1./desired_freqres; % in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_msg = '';
interpolate = 0;
sp = 1 ./ SAMPLING_FREQ;
[nd, epoch_size, nepochs] = size(data);
real_freqres = 1 ./ (sp*epoch_size);

% --- Adjust epochs if needed for frequency resolution ---
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
    [nd, epoch_size, nepochs] = size(data);
    real_freqres = 1 ./ (sp*epoch_size);
elseif round(real_freqres*100)/100 > round(desired_freqres*100)/100
    interpolate = 1;
end

% === Step 1: Channel alignment ===
if length(cnames) == length(desired_order_10_20) && all(ismember(lower(desired_order_10_20), lower(cnames)))
    % --- Strict 10-20 montage: enforce desired order ---
    % Fix common aliases (T7->T3, etc.)
    cnames = fix_aliases(cnames);
    [~, od] = ismember(lower(desired_order_10_20), lower(cnames));
    data   = data(od,:,:);
    cnames = cnames(od);
else
    % --- Extended montage: use properties.channel_params.labels ---
    if isfield(properties, 'channel_params')
        labels = properties.channel_params.labels;
    else
        labels = cnames; % fallback
    end
    labels = fix_aliases(labels);

    [~, od] = ismember(lower(labels), lower(cnames));
    data   = data(od,:,:);
    cnames = labels;

    % --- Move reference channel to last ---
    if ~strcmpi(reference, 'average')
        ref_idx = find(strcmpi(cnames, reference), 1);
        if ~isempty(ref_idx)
            cnames = [cnames(1:ref_idx-1), cnames(ref_idx+1:end), cnames(ref_idx)];
            data   = cat(1, data([1:ref_idx-1, ref_idx+1:end, ref_idx],:,:));
        end
    end
end

% === Step 2: Build output structure ===
if nepochs >= 1
    data_struct.name   = data_code;
    data_struct.srate  = SAMPLING_FREQ;
    data_struct.nchan  = length(cnames);
    data_struct.dnames = cnames;
    data_struct.ref    = reference;
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

    % === Step 3: Compute spectra & cross-spectra ===
    [Spec, data_struct.fmin, data_struct.freqres, data_struct.fmax, data_struct.CrossM, ...
        data_struct.ffteeg] = eeg_to_sp(data, sp);
    data_struct.freqrange = data_struct.fmin:data_struct.freqres:data_struct.fmax;

    % === Step 4: Interpolation if needed ===
    if interpolate
        real_frange    = data_struct.fmin:data_struct.freqres:data_struct.freqres*size(Spec,2);
        desired_frange = desired_freqres:desired_freqres:desired_freqres*size(Spec,2);
        Sp = Spec;
        for h=1:size(Spec,1)
            Spec(h,:) = pchip(real_frange, Spec(h,:), desired_frange);
        end
        Spec(:,1) = Sp(:,1);
        ii = find(Spec < 0);
        if ~isempty(ii)
            data_struct = [];
            error_msg = 'The interpolation produced negative values in the spectra';
            return
        end
        ii = find(Spec < eps);
        Spec(ii) = 1.0e-7;
        maxnormf = 49*desired_freqres;
        ii = find(desired_frange <= maxnormf);
        Spec = Spec(:, ii);
        data_struct.Spec = Spec;
        data_struct.Spec_freqrange = desired_frange(ii);
    else
        data_struct.Spec = Spec;
        data_struct.Spec_freqrange = data_struct.freqrange;
    end
else
    data_struct = [];
    error_msg = 'No epochs available';
end

end

%% === Helper to fix aliases like T7->T3, P7->T5 ===
function cnames = fix_aliases(cnames)
[ii, pos] = ismember('t7', lower(cnames)); if pos ~= 0, cnames{pos} = 'T3'; end
[ii, pos] = ismember('t8', lower(cnames)); if pos ~= 0, cnames{pos} = 'T4'; end
[ii, pos] = ismember('p7', lower(cnames)); if pos ~= 0, cnames{pos} = 'T5'; end
[ii, pos] = ismember('p8', lower(cnames)); if pos ~= 0, cnames{pos} = 'T6'; end
end

%% === Original cross-spectrum computation (unchanged) ===
function [Spec, fmin, freqres, fmax, CrossM, ffteeg] = eeg_to_sp(EEG, sp)
[nd, nt, ne] = size(EEG);
NFFT = nt;
freqres = 1 ./ (sp*NFFT);
fmin = freqres;
fmax = 1 ./ (2*sp);

ncoefs = floor(nt./2);
ffteeg = zeros(nd, ncoefs, ne);
CrossM = zeros(nd, nd, ncoefs);

for k=1:ne
    A = EEG(:,:,k);  % no Hanning window (compatibility with Cuban Norms)
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

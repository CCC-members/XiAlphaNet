function [data_struct,error_msg] = ImportEEG(properties,file_name,SubID,pat)
% IMPORTEEG
% Unified EEG import function.
% All formats are converted to a canonical PLG-like representation:
%   data ∈ R^{Nc × (lwin·fs) × Nepochs}
% before spectral processing.
%
% This guarantees identical frequency resolution and frequency range
% across file formats.

import functions.*
import plugins.*
import plugins.PLGReader.*
import plugins.HarMNqEEG.*
import functions.import.*

error_msg = '';

% -------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------
[path,~,ext] = fileparts(file_name);
ext   = lower(ext);

state        = properties.general_params.data.state;
lwin         = properties.general_params.data.lwin;   % window length in seconds (e.g. 2.56)
country      = properties.general_params.data.country;
eeg_device   = properties.general_params.data.eeg_device;
reference    = properties.general_params.data.reference;
keep_signal  = properties.general_params.data.keep_signal;

% -------------------------------------------------------------
% 1) Load data (format-specific, minimal logic only)
% -------------------------------------------------------------
switch ext

    case '.plg'
        [data,MONTAGE,~,fs,~,~,~] = read_plgwindows(file_name,state,lwin);
        data = data(1:19,:,:);

        dnames = strrep(num2cell(MONTAGE(1:19,1:3),2),'_','');
        if isequal(str2num(char(dnames)),(1:19)')
            dnames = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
                      'F7','F8','T3','T4','T5','T6','Fz','Cz','Pz'};
        end

        EEG.data     = data;
        EEG.srate    = fs;
        EEG.chanlocs = cell2struct(dnames,'labels',2);

    case '.set'
        EEG = pop_loadset(file_name);

    case '.edf'
        EEG = pop_biosig(file_name);

    case '.vhdr'
        [base_path,hdrfile,extf] = fileparts(file_name);
        EEG = pop_loadbv(base_path,[hdrfile extf]);

    case '.mat'
        S = load(file_name);   % expects variables: data, labels
        EEG = eeg_emptyset;
        EEG.data     = S.data;
        EEG.srate    = 256;
        EEG.chanlocs = cell2struct(S.labels,'labels',2);

    case '.txt'
        EEG = eeg_emptyset;
        EEG.data  = readmatrix(file_name)';
        EEG.srate = 200;

    otherwise
        error('Unsupported file format: %s',ext)
end

% -------------------------------------------------------------
% 2) Canonical PLG-style epoching (THE key step)
% -------------------------------------------------------------
EEG = standardize_epochs_PLG(EEG,lwin);

% -------------------------------------------------------------
% 3) Metadata
% -------------------------------------------------------------
dnames = {EEG.chanlocs.labels};
data   = EEG.data;
fs     = EEG.srate;

Age = '';
Sex = '';

if isfield(pat,'age'),        Age = pat.age;
elseif isfield(pat,'Age'),    Age = pat.Age; end

if isfield(pat,'sex'),        Sex = pat.sex;
elseif isfield(pat,'Sex'),    Sex = pat.Sex;
elseif isfield(pat,'Gender'), Sex = pat.Gender; end

% -------------------------------------------------------------
% 4) Spectral processing (now format-invariant)
% -------------------------------------------------------------
[data_struct,error_msg] = data_gatherer_v2( ...
    data, fs, dnames, SubID, reference, Age, Sex, ...
    country, eeg_device, keep_signal);

% -------------------------------------------------------------
% 5) Optional subject fields
% -------------------------------------------------------------
if isfield(pat,'BirthDate'),  data_struct.BirthDate  = pat.BirthDate; end
if isfield(pat,'RecordDate'), data_struct.RecordDate = pat.RecordDate; end
if isfield(pat,'Status'),     data_struct.Status     = pat.Status; end

function EEG = standardize_epochs_PLG(EEG,lwin)
% STANDARDIZE_EPOCHS_PLG
% Converts any EEG structure into PLG-style epochs:
%   data ∈ R^{Nc × (lwin·fs) × Nepochs}

data = double(EEG.data);
fs   = EEG.srate;

% Flatten if already epoched
if ndims(data) == 3
    data = reshape(data,size(data,1),[]);
end

nt = round(lwin * fs);
nw = floor(size(data,2)/nt);

if nw < 2
    error('Not enough data for %.2f s epochs.',lwin);
end

data = data(:,1:nw*nt);
data = reshape(data,size(data,1),nt,nw);

EEG.data   = data;
EEG.pnts  = nt;
EEG.trials = nw;
EEG.xmin  = 0;
EEG.xmax  = (nt-1)/fs;
EEG.times = (0:nt-1)/fs * 1000;
end


end

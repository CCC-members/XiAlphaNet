function data_struct = ImportEEG(properties,file_name,SubID,pat)

import functions.*
import plugins.*
import plugins.PLGReader.*
import plugins.HarMNqEEG.*
import functions.import.*

[path,name,ext] = fileparts(file_name);
state = properties.general_params.data.state;
lwin = properties.general_params.data.lwin;
country = properties.general_params.data.country;
eeg_device = properties.general_params.data.eeg_device;
ref = properties.general_params.data.reference;
keep_signal = properties.general_params.data.keep_signal;

switch ext
    case '.PLG'   
        [data, MONTAGE, age, SAMPLING_FREQ, epoch_size, wins, msg] = read_plgwindows(file_name, state, lwin);               
        nt                          = lwin*SAMPLING_FREQ;
        nw                          = size(data,2)./nt;
        data                        = reshape(data(1:19,1:nw*nt), 19, nt, nw);
        dnames                      = strrep(num2cell(MONTAGE(1:19, 1:3),2), '_', '');
        xx                          = str2num(char(dnames));
        if (isequal(xx, [1:19]'))
            dnames = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'F7', 'F8', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz'};
        end
        pat                         = read_plgpat(strrep(file_name, ext, '.pat')); 
        srate                       = SAMPLING_FREQ;
        Age                         = pat.Age;
        Sex                         = pat.Sex;     
    case '.edf'
        EEG                     = pop_biosig(file_name);
        EEG.setname             = SubID;
        EEG.SubID               = SubID;        
        EEG.filename            = path;
        EEG.filepath            = name;        
        % For cuban dataset
        new_labels              = replace({EEG.chanlocs.labels}','-REF','');
        [EEG.chanlocs.labels]   = new_labels{:};
        new_labels              = replace({EEG.chanlocs.labels}',' ','');
        [EEG.chanlocs.labels]   = new_labels{:};
        dnames                  = {EEG.chanlocs(1:19).labels};
        indx = ismember(dnames,{'T7','T8','P7','P8'});
        dnames(indx) = {'T3','T4','T5','T6'};
        data = EEG.data(1:length(dnames),:); 
        srate                       = EEG.srate;
        if(isfield(pat,'age' ))
            Age                         = pat.age;
        else
            Age                         = '';
        end
        Sex                         = pat.sex;
    case '.set'
         EEG                     = pop_loadset(file_name);
    case '.mat'
        EEG                     = eeg_emptyset;
        load(file_name);
        EEG.data                = data;
        
        % For Pedrito's data selection
        %         srate                   = SAMPLING_FREQ;
        %         EEG.srate               = srate;
        %         EEG.age                 = age;
        %         if(exist('labels','var'))
        %             EEG.chanlocs(length(labels)+1:end,:)    = [];
        %             new_labels                              = labels;
        %             [EEG.chanlocs.labels]                   = new_labels{:};
        %         end

        % For DEAP dataset 
        EEG.srate               = 256;  
        EEG.trials              = 1;
        EEG.nbchan              = size(data,1);
        EEG.pnts                = size(data,2);
        EEG.xmin                = 0;
        EEG.xmax                = EEG.xmin+(EEG.pnts-1)*(1/EEG.srate);
        EEG.times               = (0:EEG.pnts-1)/EEG.srate.*1000;           
        EEG.chanlocs                   = cell2struct(labels, 'labels',2);
    case '.txt'
        EEG                     = eeg_emptyset;
        [~,filename,~]          = fileparts(file_name);
        EEG.filename            = filename;
        EEG.filepath            = filepath;
        EEG.subject             = subID;
        data                    = readmatrix(file_name);
        data                    = data';
        EEG.data                = data;
        EEG.nbchan              = length(EEG.chanlocs);
        EEG.pnts                = size(data,2);
        EEG.srate               = 200;
        EEG.min                 = 0;
        EEG.max                 = EEG.xmin+(EEG.pnts-1)*(1/EEG.srate);
        EEG.times               = (0:EEG.pnts-1)/EEG.srate.*1000;
end

[data_struct, error_msg]    = data_gatherer_v2(data, srate, dnames, SubID, ref, Age, Sex,...
    country, eeg_device, keep_signal);
if(isfield(pat,'BirthDate'))
    data_struct.BirthDate       = pat.BirthDate;
end
if(isfield(pat,'RecordDate'))
    data_struct.RecordDate      = pat.RecordDate;
end
if(isfield(pat,'Status'))
    data_struct.Status          = pat.Status;
end

end


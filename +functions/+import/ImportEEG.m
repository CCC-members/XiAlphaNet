function data_struct = ImportEEG(properties,file_name)

import functions.*
import plugins.*
import plugins.PLGReader.*
import plugins.HarMNqEEG.*

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
        [nd, nt, nv]                = size(data);
        dnames                      = strrep(num2cell(MONTAGE(1:19, 1:3),2), '_', '');
        xx                          = str2num(char(dnames));
        if (isequal(xx, [1:19]'))
            dnames = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'F7', 'F8', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz'};
        end
        pat                         = read_plgpat(strrep(file_name, ext, '.pat'));    
        [data_struct, error_msg]    = data_gatherer_v2(data, SAMPLING_FREQ, dnames, pat.Name, ref, pat.Age, pat.Sex,...
            country, eeg_device, keep_signal);
        data_struct.BirthDate       = pat.BirthDate;
        data_struct.RecordDate      = pat.RecordDate;
        data_struct.Status          = pat.Status;        
    case '.set'

    case '.mat'

    case '.txt'

end


end


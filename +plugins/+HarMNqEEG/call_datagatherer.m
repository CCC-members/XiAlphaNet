%pathd = 'C:\Users\Lenovo\Desktop\lilia_raw2\all_pre_files\';
pathd = 'H:\Datos\Peter\MaLuisa\Casos-ParkinsonNeuroEpo\temp\neuroepopre\CARIDAD CABRERA RDG\';

x = dir(fullfile(pathd, '*.plg'));

data_code = {};
reference = {};
age = {};
sex = {};
country = {};
eeg_device = {};
keep_signal = 1;

addpath('H:\_VIP\Work\BVA\qeeg_tools')
ref = 'LinkedEars';
pais = 'Neuroepo';
equipo = 'MEDICID-4';
path_save = 'H:\Datos\Peter\MaLuisa\Casos-ParkinsonNeuroEpo\temp\DataGathered\neuroepopre\';
mat_sufix = '_cross_Neuroepo.mat';
state = 'A';
lwin = 2.56; %sec

cnt = 1;

for k = 1:length(x)
    fsets{k} = [pathd x(k).name];
    
    [data, MONTAGE, age, SAMPLING_FREQ, epoch_size, wins, msg] = read_plgwindows(fsets{k}, state, lwin);
    if ~isempty(ss) || isempty(data)
        disp(['Ignored: ' fsets{k}]);
        continue
    end
    disp(fsets{k})
    nt = lwin*SAMPLING_FREQ;
    nw = size(data,2)./nt;
    data = reshape(data(1:19,1:nw*nt), 19, nt, nw);
    [nd, nt, nv] = size(data);
    dnames = strrep(num2cell(MONTAGE(1:19, 1:3),2), '_', '');
    xx = str2num(char(dnames));
    if (isequal(xx, [1:19]'))
        dnames = s1020_dnames;
    end
    pat = read_plgpat(strrep(lower(fsets{k}), '.plg', '.pat'));
    names{k} = pat.Name;
    sex = pat.Sex;
    edad = age;
    [pp st ee] = fileparts(fsets{k});
    [data_struct, error_msg] = data_gatherer_v2(data, SAMPLING_FREQ, dnames, st, ref, edad, sex,...
        pais, eeg_device, keep_signal);
    if ~isempty(data_struct)
        save([path_save st mat_sufix], 'data_struct')
    end
end

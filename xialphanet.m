function xialphanet(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Automatic EEG Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies
% EGGLAB latest version
%

%%
%% Preparing WorkSpace
%%
clc;
% close all;
% restoredefaultpath;
% clearvars -except varargin;
disp('-->> Starting process');
disp('==========================================================================');

%%
%% Preparing properties
%%
addpath(genpath(fullfile(pwd,'app')));
addpath(fullfile(pwd,'dependencies'));
addpath(genpath(fullfile(pwd,'functions')));
addpath(genpath(fullfile(pwd,'guide')));
addpath(genpath(fullfile(pwd,'report')));

addpath(genpath(fullfile(pwd,'tools')));

%%
%% Init processing
%%
try
    properties = jsondecode(fileread(fullfile(pwd,'app','properties.json')));
catch EM
    fprintf(2,"\n ->> Error: The app/properties file do not have a correct format \n");
    disp("-->> Message error");
    disp(EM.message);
    disp('-->> Process stopped!!!');
    return;
end
%% Printing data information
disp(strcat("-->> Name:",properties.generals.name));
disp(strcat("-->> Version:",properties.generals.version));
disp(strcat("-->> Version date:",properties.generals.version_date));
disp("==========================================================================");

%%
%% Starting mode
%%
setGlobalGuimode(true);
for i=1:length(varargin)
    if(isequal(varargin{i},'nogui'))
        setGlobalGuimode(false);
    end
end
if(getGlobalGuimode())
    XiAlphaNet
else
    properties = get_properties();
    [status, properties] = check_properties(properties);
    if(status)
        process_interface(properties);
    end

    %% Finishing process
    disp('==========================================================================');
    disp('--------------------------Process finished.-------------------------------');
    disp('==========================================================================');
    disp(strcat("-->> ",properties.generals.name));
    disp('==========================================================================');
    close all;
    clear all;

end

end

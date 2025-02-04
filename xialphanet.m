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
disp(strcat("-->> Name:",properties.app.name));
disp(strcat("-->> Version:",properties.app.version));
disp(strcat("-->> Version date:",properties.app.version_date));
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
    %checked = check_properties()
    processing_interface(properties)

    %% Finishing process
    disp('--------------------------Process finished.-------------------------------');
    disp('==========================================================================');
    disp(properties.generals.name);
    disp('==========================================================================');
    close all;
    clear all;

end

end

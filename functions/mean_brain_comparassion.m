% This program extracts the activation and peak frequency distribution for
% Xi and Alpha processes in the brain and creates plots to visually compare
% the differences between Control and Pathological groups.
% The first row in the plot correspond to the Control group, the
% second to the pathological group while the last row correspond to the
% Absolute  difference between groups
clc;clear all;

modelParametersPath = 'Data\Model_Parameters';
dataAge =  'Data\Age';

% Subfolders within the main folder
subFolders = {'Control'};%, 'Pathological'};

% 
load('Data\Model_Parameters\parameters.mat');
Ne = parameters.Dimensions.Ne;
Nr = parameters.Dimensions.Nr;
Nv = parameters.Dimensions.Nv;
Nw = parameters.Dimensions.Nw;
x0 = zeros(7*Nv+1,1);
for k = 1:length(subFolders)
    % Full path to the subfolder

    folderPath_X = fullfile(modelParametersPath, subFolders{k});
    folderPath_Age = fullfile(dataAge, subFolders{k});
    
    % Get a list of all .mat files in the subfolder
    matFiles_X = dir(fullfile(folderPath_X, '*.mat'));
    matFiles_Age = dir(fullfile(folderPath_Age, '*.mat'));
   
    %% Loop through each .mat file in the subfolder
    for j = 1:length(matFiles_X)
        %
        filePath_X = fullfile(folderPath_X, matFiles_X(j).name);
        % Load the .mat file
        x = load(filePath_X);
        x = x.x;
        %[x_opt] = np_ref_solution(x);
        x = x.Solution;
        x0(:,k) = x0(:,k) + x;
    end
    x0(:,k)=x0(:,k)/length(matFiles_X);
end

[e1,a1,s1] = x2v(x0(:,1));
% a2(:,4)=exp(a2(:,4));
% a1(:,4)=exp(a1(:,4));
% 
% colorlimit_e_max = max(max(e1,e2));
% colorlimit_a_max = max(max(a1,a2));
% colorlimit_e_min = min(min(e1,e2));
% colorlimit_a_min = min(min(a1,a2));
% 
% subplot(3,3,1)
% Ja1= a1(:,1);
% colorlimit_alpha_activation = [colorlimit_a_min(1),colorlimit_a_max(1)];
% esi_plot(gca,Ja1,colorlimit_alpha_activation);
% 
% 
% subplot(3,3,2)
% Je1= e1(:,1);
% colorlimit_xi_activation = [colorlimit_e_min(1),colorlimit_e_max(1)];
% esi_plot(gca,Je1,colorlimit_xi_activation);
% 
% 
% subplot(3,3,3)
% Ja14= a1(:,4);
% colorlimit_alpha_freq = [colorlimit_a_min(4),colorlimit_a_max(4)];
% esi_plot(gca,Ja14,colorlimit_alpha_freq);
% 
% 
% subplot(3,3,4)
% Ja2= a2(:,1);
% colorlimit_alpha_activation = [colorlimit_a_min(1),colorlimit_a_max(1)];
% esi_plot(gca,Ja2,colorlimit_alpha_activation);
% 
% 
% subplot(3,3,5)
% Je2= e2(:,1);
% colorlimit_xi_activation = [colorlimit_e_min(1),colorlimit_e_max(1)];
% esi_plot(gca,Je2,colorlimit_xi_activation);
% 
% 
% subplot(3,3,6)
% Ja24= a2(:,4);
% colorlimit_alpha_freq = [colorlimit_a_min(4),colorlimit_a_max(4)];
% esi_plot(gca,Ja24,colorlimit_alpha_freq);
% 
% 
% a2(:,1) = threshold_to_zero(a2(:,1));
% a1(:,1) = threshold_to_zero(a1(:,1));
% 
% a2(:,4) = threshold_to_zero(a2(:,4));
% a1(:,4) = threshold_to_zero(a1(:,4));
% 
% e2(:,1) = threshold_to_zero(e2(:,1));
% e1(:,1) = threshold_to_zero(e1(:,1));
% 
% subplot(3,3,7)
% Ja2= abs(a2(:,1)-a1(:,1));
% colorlimit_alpha_activation = [colorlimit_a_min(1),colorlimit_a_max(1)];
% esi_plot(gca,Ja2);
% 
% 
% subplot(3,3,8)
% Je2= abs(e2(:,1)-e1(:,1));
% colorlimit_xi_activation = [colorlimit_e_min(1),colorlimit_e_max(1)];
% esi_plot(gca,Je2);
% 
% 
% subplot(3,3,9)
% Ja24= abs(a2(:,4)-a1(:,4));
% colorlimit_alpha_freq = [colorlimit_a_min(4),colorlimit_a_max(4)];
% esi_plot(gca,Ja24);
% 
% 
% colormap(turbo)




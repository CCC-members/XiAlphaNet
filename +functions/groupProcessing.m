function XIALPHANET = groupProcessing(properties,XIALPHANET,varargin)

import functions.auxx.ModelVectorization.*
age_range = properties.model_params.group.age_range;

age_min = 0;%age_range(1);
age_max = 100;%age_range(2);
parameters = load(fullfile(XIALPHANET.Location,XIALPHANET.Structural.parameters));
ages = [];
All_Data = {}; 
index = 1;
for i=1:length(XIALPHANET.Participants)
    participant = XIALPHANET.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed') && age_min <= participant_age && participant_age <= age_max)
       ages = [ages,participant_age];
       All_Data{2,index} =  participant_age;
       Part_Info = jsondecode(fileread(fullfile(XIALPHANET.Location,participant.SubID,participant.FileInfo)));
       alpha_process = load(fullfile(XIALPHANET.Location,participant.SubID,Part_Info.Alpha_estimate));
       a(:,1) = alpha_process.Power;
       a(:,2) = alpha_process.Width;
       a(:,3) = alpha_process.Exponent;
       a(:,4) = alpha_process.PAF;
       xi_process = load(fullfile(XIALPHANET.Location,participant.SubID,Part_Info.Xi_estimate));
       e(:,1) = xi_process.Power;
       e(:,2) = xi_process.Width;
       e(:,3) = xi_process.Exponent;
       s2 = 1;
       x = v2x(e,a,s2);
       All_Data{1,index} = x;
       index = index +1;
       % if (isequal(Process,'Alpha'))
       %     spec_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
       % elseif (isequal(Process,'Xi'))
       %     spec_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
       % end
       % if isequal(variable,'Power')
       %      All_Data{1,i} = spec_process.Power;
       % elseif (isequal(variable,'Width'))
       %      All_Data{1,i} = spec_process.Width;
       % elseif (isequal(variable,'Exponent'))
       %      All_Data{1,i} = spec_process.Exponent;
       % elseif (isequal(variable,'PAF')) && (isequal(Process,'Alpha'))
       %      All_Data{1,i} = spec_process.PAF;
       % end
    end
end


%% Simple Counting 
% Initialize storage for Peak Alpha Frequency (a(:,4)), Amplitude of the Alpha (a(:,1)), and Amplitude of Xi
PAF_all = [];
AlphaAmp_all = [];
XiAmp_all = [];

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    PAF_all(:,j) = a(:,4);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = a(:,1);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = e(:,1);                                % Store Amplitude of the Xi
end


% Define threshold for each measure use set_threshold_em for recalculate this value
threshold_PAF = 7;             % Example threshold for PAF
% threshold_AlphaAmp = 0.19;   % Example threshold for Alpha Amplitude
% threshold_XiAmp = 0.19;      % Example threshold for Xi Amplitude

% Convert each data.x to [e, a, s2] and store PAF (a(:,4)), Alpha Amplitude (a(:,1)), and Xi Amplitude
for j = 1:length(All_Data(1,:))
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold_AlphaAmp   =  prctile(AlphaAmp_all(:,j),75);
    PAF_all(:,j) = PAF_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);            % Store Peak Alpha Frequency
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j).*(AlphaAmp_all(:,j)>threshold_AlphaAmp);       % Store Amplitude of the Alpha
    XiAmp_all(:,j) = XiAmp_all(:,j);         % Store Amplitude of the Xi
end


% Apply threshold and create binary masks for each measure
for j = 1:length(All_Data(1,:))
    fprintf('Processing subject %d\n', j);
    %[e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold_AlphaAmp   =  prctile(AlphaAmp_all(:,j),75);
    PAF_all(:,j) = PAF_all(:,j) > threshold_PAF;             % Apply threshold for PAF
    AlphaAmp_all(:,j) = AlphaAmp_all(:,j) > threshold_AlphaAmp;   % Apply threshold for Alpha Amplitude
    XiAmp_all(:,j) =  XiAmp_all(:,j) > threshold_AlphaAmp;         % Apply threshold for Xi Amplitude
end


% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages

% Define age intervals
age_intervals = linspace(age_min, age_max, floor((age_max-age_min)/20)+1);  % 5 intervals, 6 edges

% Initialize storage for the averaged values
PAF_avg_intervals = zeros(size(PAF_all, 1), length(age_intervals)-1);
AlphaAmp_avg_intervals = zeros(size(AlphaAmp_all, 1), length(age_intervals)-1);
XiAmp_avg_intervals = zeros(size(XiAmp_all, 1), length(age_intervals)-1);

% Average values for each age interval
for i = 1:length(age_intervals)-1
    age_mask = (ages >= age_intervals(i)) & (ages < age_intervals(i+1));
    
    % For PAF
    f_va_PAF = sum(PAF_all(:, age_mask), 2);
    nf_va_PAF = sum(age_mask);  % Number of subjects in this age group
    PAF_avg_intervals(:, i) = f_va_PAF / (1+nf_va_PAF);
    
    % For Alpha Amplitude
    f_va_AlphaAmp = sum(AlphaAmp_all(:, age_mask), 2);
    nf_va_AlphaAmp = sum(age_mask);  % Number of subjects in this age group
    AlphaAmp_avg_intervals(:, i) = f_va_AlphaAmp / (1+nf_va_AlphaAmp);
    
    % For Xi Amplitude
    f_va_XiAmp = sum(XiAmp_all(:, age_mask), 2);
    nf_va_XiAmp = sum(age_mask);  % Number of subjects in this age group
    XiAmp_avg_intervals(:, i) = f_va_XiAmp / (1+nf_va_XiAmp);
end

[XIALPHANET] = xan_save(properties,[],'groups',XIALPHANET,data1,data2,data3);


end


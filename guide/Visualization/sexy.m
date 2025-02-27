clc;
clear all

age_min = 0;%age_range(1);
age_max = 20;%age_range(2);
dataset = jsondecode(fileread('/home/ronaldo/Documents/dev/Data/Results/XIALPHANET.json'));
parameters = load('/home/ronaldo/Documents/dev/Data/Results/structural/parameters.mat');
ages = [];
All_Data = {}; 
index = 1;
for i=1:length(dataset.Participants)
    participant = dataset.Participants(i);
    participant_age = participant.Age;
    if(isequal(participant.Status,'Completed')) && age_min <= participant_age && participant_age <= age_max
       ages = [ages,participant_age];
       All_Data{2,index} =  participant_age;
       Part_Info = jsondecode(fileread(fullfile(dataset.Location,participant.SubID,participant.FileInfo)));
       alpha_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Alpha_estimate));
       a(:,1) = alpha_process.Power;
       a(:,2) = alpha_process.Width;
       a(:,3) = alpha_process.Exponent;
       a(:,4) = alpha_process.PAF;
       xi_process = load(fullfile(dataset.Location,participant.SubID,Part_Info.Xi_estimate));
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


%%
% Initialize storage for a(:,4) for all subjects
a_all = [];

% Convert each data.x to [e, a, s2] and store a(:,4)
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    a_all(:,j) = a(:,1);  % Store the fourth column of a in a_all
end
%
for j = 1:length(All_Data(1,:))
    [e, a, s2] = x2v(All_Data{1,j});  % Assuming x2v returns [e, a, s2]
    threshold = prctile(a_all(:,j),80);
    a_all(:,j) = a(:,1)>threshold;  % Store the fourth column of a in a_all
end

% Get unique ages
ages = cell2mat(All_Data(2,:));  % Get the ages
unique_ages = unique(ages);

% Initialize storage for the averaged values
a_avg = zeros(size(a_all, 1), length(unique_ages));

% Average a(:,4) values for the same age
for i = 1:length(unique_ages)
    age_mask = (ages == unique_ages(i));
    % Sum for each voxel the number of occurrences across the subjects with the same age
    f_va = sum(a_all(:, age_mask), 2);
    % Normalize by the number of subjects in this age group
    nf_va = sum(age_mask);  % Number of subjects in this age group
    a_avg(:, i) = f_va / nf_va;
end

% Apply a moving average to smooth the transitions
windowSize = 5;  % Define the size of the moving average window (adjust as needed)
a_smooth = movmean(a_avg, windowSize, 2);

% Determine global color limits across all ages
global_min = min(a_smooth(:));
global_max = max(a_smooth(:));
colorLimits = [0, global_max];

% Initialize Video Writer
videoFilename = 'J_age_video.avi';
v = VideoWriter(videoFilename, 'Uncompressed AVI');
v.FrameRate = 10;
open(v);

% Create the video by plotting J(age) = a_smooth(:,i) for each unique age
figure;
ax = axes;  % Create an axis for plotting

for i = 1:length(unique_ages)
    i
    J_age = a_smooth(:,i);  % Get the smoothed a(:,4) at this age
    
    % Use the esi_plot function for plotting with fixed color limits
    esi_plot(ax, J_age, colorLimits);
    %colormap("hot");
    title(ax, sprintf('Age: %.2f', unique_ages(i)));
    
    % Capture the frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Finalize the Video
close(v);

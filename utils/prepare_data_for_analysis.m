%% ========================================================================
% EEG DATASET ORGANIZER & METADATA EXTRACTOR
%
% PURPOSE:
%   - Organize EEG dataset by session & subject.
%   - Copy only files matching a given task (e.g., "eyesclosed").
%   - Rename copied files into short form (last underscore part).
%   - Extract demographic metadata (Age, Sex) if available in .mat file.
%   - Save Participants.json per session.

% INPUTS (edit these before running):
%   inputRoot  = path to original dataset
%   outputRoot = path to save organized dataset
%   task_class = string pattern of task (e.g., "eyesclosed")
%
% DEPENDENCIES:
%   - Requires saveJSON (from your tools package).
%

%% ========================================================================
 
%% --- USER CONFIGURATION ---

inputRoot  = 'E:\A test-retest resting and cognitive state EEG dataset';
outputRoot = 'E:\OrganizedDataset';
task_class = 'eyesclosed';  % Example: "eyesclosed", "oddball", etc.

% Ensure output directory exists

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end
 
%% --- STEP 1: DETECT SUBJECTS ---

subjects = dir(inputRoot);
subjects(ismember({subjects.name}, {'.','..'})) = [];  % remove system entries
subjects(~[subjects.isdir]) = [];                      % keep only directories
fprintf('Found %d subjects in %s\n', numel(subjects), inputRoot);
 
%% --- STEP 2: LOOP OVER SESSIONS ---

% Collect all unique session names across subjects

allSessions = {};
for iSub = 1:numel(subjects)
    subjPath = fullfile(inputRoot, subjects(iSub).name);
    sesDirs = dir(subjPath);
    sesDirs(ismember({sesDirs.name}, {'.','..'})) = [];
    sesDirs(~[sesDirs.isdir]) = [];
    allSessions = [allSessions, {sesDirs.name}];
end

allSessions = unique(allSessions);
 
fprintf('Detected sessions: %s\n', strjoin(allSessions, ', '));
 
%% --- STEP 3: PROCESS EACH SESSION ---

for iSes = 1:numel(allSessions)

    sesName = allSessions{iSes};
    sesOutputDir = fullfile(outputRoot, sesName);
    if ~exist(sesOutputDir, 'dir')
        mkdir(sesOutputDir);
    end 
    fprintf('\n-->> Processing session: %s\n', sesName);
    % Initialize metadata container for this session
    Participants = struct([]);
    % --- Loop over subjects ---
    for iSub = 1:numel(subjects)
        SubID = subjects(iSub).name;
        subjPath = fullfile(inputRoot, SubID, sesName, 'eeg');
        if ~exist(subjPath, 'dir')
            fprintf('   Skipping %s - no eeg folder for %s.\n', SubID, sesName);
            continue;
        end
        % Create subject output folder
        subjOutputDir = fullfile(sesOutputDir, SubID);
        if ~exist(subjOutputDir, 'dir')
            mkdir(subjOutputDir);
        end
        % --- Copy task files ---
        files = dir(fullfile(subjPath, '*'));
        copiedAny = false;
        for f = 1:length(files)
            if ~files(f).isdir && contains(files(f).name, ['_task-' task_class '_'])
                [~, name, ext] = fileparts(files(f).name);
                parts = strsplit(name, '_');
                shortName = [parts{end}, ext];
                srcFile = fullfile(subjPath, files(f).name);
                dstFile = fullfile(subjOutputDir, shortName);
                copyfile(srcFile, dstFile);
                fprintf('   Copied: %s → %s\n', files(f).name, shortName);
                copiedAny = true;
            end
        end
 
        if ~copiedAny
            fprintf('   No task "%s" files found for %s in %s.\n',task_class, SubID, sesName);
        end
 
        % --- Extract metadata if available ---

        matFile = fullfile(inputRoot, SubID, [SubID '.mat']);
        metaEntry = struct('SubID', SubID, 'Age', '', 'Sex', '');
        if exist(matFile, 'file')
            try
                loadedData = load(matFile);
                if isfield(loadedData, 'data_struct')
                    metaEntry.Age = loadedData.data_struct.age;
                    metaEntry.Sex = loadedData.data_struct.sex;
                end
                fprintf('   Metadata loaded for %s\n', SubID);
            catch
                fprintf('  Failed to load metadata for %s\n', SubID);
            end
        end
 
        % Append entry
        if isempty(Participants)
            Participants = metaEntry;
        else
            Participants(end+1) = metaEntry; 
        end
    end
 
    % --- Save JSON for this session ---

    import tools.*
    saveJSON(Participants, fullfile(sesOutputDir, 'Participants.json'));
    fprintf('   Participants.json saved in %s\n', sesOutputDir);
end
 
fprintf('\nAll sessions processed successfully.\n');

 
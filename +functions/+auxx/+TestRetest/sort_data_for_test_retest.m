% Configuration
inputRoot = 'E:\A test-retest resting and cognitive state EEG dataset';
outputRoot = 'E:\Test-retest_Eyesclosed';  % ️ CHANGE THIS TO YOUR DESIRED OUTPUT PATH

% Create output root if not exists
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

% Subjects: sub-01 to sub-59
subjects = cell(1, 59);
for i = 1:59
    subjects{i} = sprintf('sub-%02d', i);
end

% Sessions
sessions = {'ses-session1', 'ses-session2', 'ses-session3'};

% Pattern to match
ecPattern = '_task-eyeclosed_';

% Loop over sessions first
for s = 1:length(sessions)
    sessionName = sessions{s};
    sessionOutputDir = fullfile(outputRoot, ['session', num2str(s)]);
    
    % Create session output folder
    if ~exist(sessionOutputDir, 'dir')
        mkdir(sessionOutputDir);
    end
    
    fprintf('Processing %s...\n', sessionName);
    
    % Loop over subjects
    for i = 1:length(subjects)
        subj = subjects{i};
        subjInputPath = fullfile(inputRoot, subj, sessionName, 'eeg');
        
        % Skip if subject/session/eeg folder doesn't exist
        if ~exist(subjInputPath, 'dir')
            fprintf('  Skipping %s: folder not found.\n', subj);
            continue;
        end
        
        % Output path for this subject in this session
        subjOutputPath = fullfile(sessionOutputDir, subj, 'eyes_closed_data');
        if ~exist(subjOutputPath, 'dir')
            mkdir(subjOutputPath);
        end
        
        % Get list of files in eeg folder
        files = dir(fullfile(subjInputPath, '*'));
        
        % Filter and copy EC files
        copiedAny = false;
        for f = 1:length(files)
            if ~files(f).isdir && contains(files(f).name, ecPattern)
                srcFile = fullfile(subjInputPath, files(f).name);
                dstFile = fullfile(subjOutputPath, files(f).name);
                copyfile(srcFile, dstFile);
                fprintf('  Copied: %s → %s\n', files(f).name, subj);
                copiedAny = true;
            end
        end
        
        if ~copiedAny
            fprintf('  No EC files found for %s in %s.\n', subj, sessionName);
        end
    end
end

fprintf(' Done! Eyes Closed data organized by session.\n');
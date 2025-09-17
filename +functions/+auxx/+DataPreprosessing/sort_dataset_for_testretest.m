% Configuration
inputRoot = 'E:\A test-retest resting and cognitive state EEG dataset'; 
outputRoot = 'E:\Test-retest_Eyesclosed';  %  CHANGE IF NEEDED

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

%  Pattern to match — note: "eyesclosed", not "eyeclosed"
ecPattern = '_task-eyesclosed_';

% Loop over sessions first
for s = 1:length(sessions)
    sessionName = sessions{s};
    sessionOutputDir = fullfile(outputRoot, ['session', num2str(s)]);
    
    % Create session output folder
    if ~exist(sessionOutputDir, 'dir')
        mkdir(sessionOutputDir);
    end
    
    fprintf(' Processing %s...\n', sessionName);
    
    % Loop over subjects
    for i = 1:length(subjects)
        subj = subjects{i};
        subjInputPath = fullfile(inputRoot, subj, sessionName, 'eeg');
        
        % Skip if subject/session/eeg folder doesn't exist
        if ~exist(subjInputPath, 'dir')
            fprintf('   Skipping %s: folder not found.\n', subj);
            continue;
        end
        
        % Output path for this subject in this session (directly under sub-XX, no subfolder)
        subjOutputPath = fullfile(sessionOutputDir, subj);
        if ~exist(subjOutputPath, 'dir')
            mkdir(subjOutputPath);
        end
        
        % Get list of files in eeg folder
        files = dir(fullfile(subjInputPath, '*'));
        
        % Filter and copy EC files, RENAMING them to short form
        copiedAny = false;
        for f = 1:length(files)
            if ~files(f).isdir && contains(files(f).name, ecPattern)
                
                % Extract short filename: remove prefix up to last underscore before extension
                [~, name, ext] = fileparts(files(f).name);
                
                % Example: "sub-02_ses-session1_task-eyesclosed_channels" → we want "channels"
                % Split by underscores and take last part
                parts = strsplit(name, '_');
                shortName = [parts{end}, ext];  % e.g., "channels.tsv"
                
                % Build source and destination paths
                srcFile = fullfile(subjInputPath, files(f).name);
                dstFile = fullfile(subjOutputPath, shortName);
                
                % Copy file
                copyfile(srcFile, dstFile);
                fprintf('  Copied: %s → %s\n', files(f).name, shortName);
                copiedAny = true;
            end
        end
        
        if ~copiedAny
            fprintf('  No eyesclosed files found for %s in %s.\n', subj, sessionName);
        end
    end
end

fprintf('Eyes closed data organized by session with simplified filenames.\n');
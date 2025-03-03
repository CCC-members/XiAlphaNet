function [] = extractFile(zipFilename, outputDir, outputFile)
% extractFile

% Obtain the entry's output names
outputName = fullfile(outputDir, outputFile);

% Create a stream copier to copy files.
streamCopier = ...
    com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;

% Create a Java zipFile object and obtain the entries.
try
    % Create a Java file of the Zip filename.
    zipJavaFile = java.io.File(zipFilename);

    % Create a java ZipFile and validate it.
    zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);

    % Get entry
    entry = zipFile.getEntry(outputFile);

catch exception
    if ~isempty(zipFile)
        zipFile.close;
    end
    delete(cleanUpUrl);
    error(message('MATLAB:unzip:unvalidZipFile', zipFilename));
end

% Create the Java File output object using the entry's name.
file = java.io.File(outputName);

% If the parent directory of the entry name does not exist, then create it.
parentDir = char(file.getParent.toString);
if ~exist(parentDir, 'dir')
    mkdir(parentDir)
end

% Create an output stream
try
    fileOutputStream = java.io.FileOutputStream(file);
catch exception
    overwriteExistingFile = file.isFile && ~file.canWrite;
    if overwriteExistingFile
        warning(message('MATLAB:extractArchive:UnableToOverwrite', outputName));
    else
        warning(message('MATLAB:extractArchive:UnableToCreate', outputName));
    end
    return
end

% Create an input stream from the API
fileInputStream = zipFile.getInputStream(entry);

% Extract the entry via the output stream.
streamCopier.copyStream(fileInputStream, fileOutputStream);

% Close the output stream.
fileOutputStream.close;

end


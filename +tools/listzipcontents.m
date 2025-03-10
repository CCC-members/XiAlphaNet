function filelist = listzipcontents(zipFilename)
%LISTZIPCONTENTS Lists the contents of zip file.
%
%   LISTZIPCONTENTS(ZIPFILENAME) lists the archived contents of ZIPFILENAME
%
%   ZIPFILENAME is a string specifying the name of the zip file. 
%   ZIPFILENAME can include the directory name; otherwise, the file must be
%   in the current directory or in a directory on the MATLAB path.
%
%   FILENAMES = LISTZIPCONTENTS(ZIPFILENAME) lists the zip file contents
%   and returns the file names into the string cell array FILENAMES.
%
%   Unsupported zip files
%   ---------------------
%   LISTZIPCONTENTS does not support password-protected or encrypted zip 
%   archives.
%
%   Examples
%   --------
%   % List the contents of demos.zip 
%   listzipcontents('demos.zip')
%
%   B.C. Hamans, UMC St Radboud, 2011 (B.C.Hamans@rad.umcn.nl)
%
%   See also FILEATTRIB, GZIP, GUNZIP, TAR, UNTAR, UNZIP, ZIP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error(nargchk(1,1,nargin,'struct'));
% Create a Java file of the ZIP filename.
zipJavaFile  = java.io.File(zipFilename);
% Create a Java ZipFile and validate it.
zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);
% Extract the entries from the ZipFile.
entries = zipFile.getEntries;
% Initialize the file list.
filelist={};
% Loop through the entries and add to the file list.
while entries.hasMoreElements
    filelist = cat(1,filelist,char(entries.nextElement));
end
zipFile.close;
end

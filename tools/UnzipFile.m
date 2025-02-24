function unzipped = UnzipFile(folder,filename)
unzipped = false;
exampleFiles = unzip(fullfile(folder,filename),folder);
pause(1);
delete(fullfile(folder,filename));
unzipped = true;
end


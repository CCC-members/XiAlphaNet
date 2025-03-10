function downloaded = DownloadFile(url,folder,filename)
downloaded = false;
try
    disp("Downloading file");
    matlab.net.http.HTTPOptions.VerifyServerName = false;
    options = weboptions('Timeout',Inf,'RequestMethod','auto');
    downladed_file = websave(fullfile(folder,filename),url,options);
    setGlobalDownloaded(true);
    downloaded = true;
catch
    downloaded = false;
end
end


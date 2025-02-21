function DownloadFile(url,folder,filename)
disp("Downloading file");
matlab.net.http.HTTPOptions.VerifyServerName = false;
options = weboptions('Timeout',Inf,'RequestMethod','auto');
downladed_file = websave(fullfile(folder,filename),url,options);
setGlobalDownloaded(true);
end


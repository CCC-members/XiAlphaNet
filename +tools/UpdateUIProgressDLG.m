function downloaded = UpdateUIProgressDLG(dlg,type,file,fullSize)
pause(3);
getGlobalDownloaded()
while ~getGlobalDownloaded()
    try
        s = dir(file);
        size = sizeconvert(s.bytes,1);
        dlg.Value = size/fullSize;
        disp(strcat("Process value: ", num2str(dlg.Value)));
        drawnow;
    catch Ex
        disp("Error");
        disp(Ex.message);
        pause(3);
    end
end
dlg.Value = 1;
drawnow;
downloaded = true;
end


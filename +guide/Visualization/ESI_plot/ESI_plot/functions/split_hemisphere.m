function [Cortex,iHideVert]  = split_hemisphere(Cortex, hemis)
[rH, lH] = tess_hemisplit(Cortex);

if(isequal(hemis,'left'))
    iHideVert = rH;
end
if(isequal(hemis,'right'))
    iHideVert = lH;
end
if(isempty(hemis))
    iHideVert = [];
end
Cortex.Vertices(iHideVert,:)    = [];
Cortex.VertConn(iHideVert,:)    = [];
Cortex.VertConn(:,iHideVert)    = [];
Cortex.SulciMap(iHideVert,:)    = [];
Cortex.VertNormals(iHideVert,:) = [];
Cortex.Curvature(iHideVert,:)   = [];
isHideFaces                     = any(ismember(Cortex.Faces, iHideVert), 2);
Cortex.Faces(isHideFaces,:)     = [];
if(~isempty(hemis))
    Cortex.Faces                = Cortex.Faces(:,:) - length(iHideVert);
end
sAtlas = Cortex.Atlas(Cortex.iAtlas);
if isempty(sAtlas)
    return;
end
sScouts = sAtlas.Scouts;
if isempty(sScouts)
    return;
end
sScouts(ismember({sScouts.Label},{'Thalamus L','Thalamus R','Unknown L','Unknown R','VentralDC L L','VentralDC L R','VentralDC R',...
    'Ventricle inf-lat L','Ventricle inf-lat R','Ventricle lat L L','Ventricle lat L R','Ventricle lat R',...
    'White L','White R'})) = [];
if(isequal(hemis,'left'))
    sScouts(ismember({sScouts.Region},{'RU','UU'})) = [];
end
if(isequal(hemis,'right'))
    sScouts(ismember({sScouts.Region},{'LU','UU'})) = [];
end
for i=1:length(sScouts)
    sScouts(i).Vertices(ismember(sScouts(i).Vertices,iHideVert)) = [];
    if(~isempty(hemis))
        sScouts(i).Vertices = sScouts(i).Vertices - length(iHideVert);
    end
end
Cortex.Atlas(Cortex.iAtlas).Scouts = sScouts;

end


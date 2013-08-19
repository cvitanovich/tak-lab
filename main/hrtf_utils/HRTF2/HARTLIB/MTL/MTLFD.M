function index = mtlfd(dir,azimuth, elevation)

% function index = mtlfd(dir,azimuth, elevation)

[m,n_dir]=size(dir);

if (length(azimuth) ~=1 ) index=[]; return; end;
if (length(elevation) ~=1) index=[]; return; end;

ele_index=find(dir(1,:)==[ones(1,n_dir)*elevation]);
azi_index=find(dir(2,ele_index)==[ones(1,length(ele_index))*azimuth]);

index=ele_index(azi_index);
return;

end;

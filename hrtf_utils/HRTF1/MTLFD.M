function index = mtlfd(dir,azimuth, elevation)

% function index = mtlfd(dir,azimuth, elevation)

[m,n] = size(dir);

ele_index=find(dir(1,:)==[ones(1,n)*elevation]);
azi_index=find(dir(2,ele_index)==[ones(1,length(ele_index))*azimuth]);
index=ele_index(azi_index);
return;

end;

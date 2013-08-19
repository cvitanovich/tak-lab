function wrmap(MAP, filename)

% writes MAP in filename 
% PC format, float , elevation azimmuth ! elevation azimuth ! .....


[fid,message] = fopen(filename,'w','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

[m,n] = size(MAP);

for i=1:n
 fwrite(fid, [MAP(1,i)  MAP(2,i)], 'float32');
end

fclose(fid);

return

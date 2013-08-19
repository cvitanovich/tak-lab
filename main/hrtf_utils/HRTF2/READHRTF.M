function [spec1, spec2, az, elev, sr, sf, ef] = readhrtf(filename);
% READHRTF  function to read an hrtf file and get the data
% [spec1, spec2, az, elev, sr, sf, ef] = readhrtf(filename)
% spec1 = complex spectrum of first data set
% spec2 = complex spectrum of second data set
% az = azimuth
% elev = elevation
% sr = sampling rate
% sf = starting freq
% ef = ending freq
% Code reads some other things that may be useful in the future.

fid = fopen(filename,'r');
if ~(fid>0) error('unable to open hrtf file'); end;

fseek(fid,128,'bof');
nlines = fread(fid,1,'short');

fseek(fid,132,'bof');
sr = (fread(fid,1,'float'));

fseek(fid,136,'bof');
begin_frq = (fread(fid,1,'float'));

fseek(fid,140,'bof');
end_frq = (fread(fid,1,'float'));

fseek(fid,222,'bof');
comment = setstr(fread(fid,200,'char')');
disp(['Comment:  ',comment]);

fseek(fid,512,'bof');
elev = (fread(fid,1,'float'))

fseek(fid,516,'bof');
az = (fread(fid,1,'float'))

fseek(fid,520,'bof');
n_avg = fread(fid,1,'short');

fseek(fid, 522, 'bof');
attenuation = fread(fid,1,'float');

if (strcmp(filename(10:12),'802'))
	sr = 24000;
	sf = 1000;
	ef = 11000;
	n_avg = 100;
	attenuation = 40;
end;
lines = 1024;
start_data = 3*256;

% read the spectrum data
fseek(fid, start_data,'bof');
for i=1:1024
	spec1(i) = fread(fid,1,'float');
	spec1(i) = spec1(i) + sqrt(-1)*fread(fid,1,'float');
end;

for i=1:1024
	spec2(i) = fread(fid,1,'float');
	spec2(i) = spec2(i) + sqrt(-1)*fread(fid,1,'float');
end;

fclose(fid);
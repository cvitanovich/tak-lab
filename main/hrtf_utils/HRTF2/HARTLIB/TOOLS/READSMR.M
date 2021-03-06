function T=readsmr(fname)

% READSMR(fname) reads smr files generated by sigmaplot.
% 
% Klaus Hartung
% 7.5.96
% Ruhr-Universitaet Bochum
% Lehrstuhl fuer allgemeine Elektrotechnik und Akustik

[fid,message] = fopen(fname);
if (fid < 0) error(message); end

counter=1;

while 1
	line=fgetl(fid);
	if (line(1)==(-1)) break; end
	if (~isempty(line)) 
		line=fliplr(deblank(fliplr(line)));
		line=strrep(line,',','.');
		T(counter)=str2num(line);
		counter=counter+1;
	end
end 
stat=fclose(fid)

%T=T.';

if (stat<0) error('readsmr: error at fclose'); end;

function [Lnoi, Rnoi] = read_mls_stim(el,az,stim_dir)

%function [Lnoi, Rnoi] = read_mls_stim(el,az)
% reads stimuli generated by mls_style programs

if nargin < 1 | isempty(el)
el = input('enter elevation : ');
end
if nargin < 2 | isempty(az)
az = input('enter azimuth : ');
end
if nargin < 3 | isempty(stim_dir)
  stim_dir = 'e:\kip\temp\ILA2\';
end
  
% load stimuli
fname = ['stimbuf_left_' num2str(el) '_' num2str(az)];
fid = fopen([stim_dir fname],'r');
[Lnoi, count] = fread(fid, inf,'int16');
fclose(fid);

fname = ['stimbuf_right_' num2str(el) '_' num2str(az)];
fid = fopen([stim_dir fname],'r');
[Rnoi, count] = fread(fid, inf,'int16');
fclose(fid);

function delay=mtlrdel(filename,n_channels)

% delay=mtlrdel(filename,n_channels)
% 
% MTLRDEL reads delay from MTLSIG catalogue
%
%         delay: array with delay in lines for each channel
%         filename: string with filename
% 
% optional parameters:
%
%         n_channels: number of channels to read

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 21.4.1995


[file_type,info_blocks,n_channel,n_lines, sampling_rate,...
 start,stop,n_directions,comment1, comment2] = mtlrh(filename);

if (nargin < 2) n_channels=n_channel; end;
if (n_channels > n_channel) error('ERROR: mtlrdel2,  n_channels greater than available channels' );
delay=-1; return; end;

[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
   disp(message);
   delay = -1;
   return;
end;

status=fseek(fid,512+n_channels*4,'bof');
if (status<0) error('ERROR: mtlrdel2,  fseek to delay not succesful' ); delay=-1; return; end;

delay=fread(fid,n_channels,'float32');

status= fclose(fid);
if (status<0) error('ERROR: mtlrdel2,  fclose not succesful' ); delay=-1; return; end;

return;


end;
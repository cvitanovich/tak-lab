function channel=mtlrch(filename, index);

% MTLRCH reads one channel from a file (type 2 or 3)
%
% [channel]=mtlrch(filename, index);
%
%         filename: string with filename
%         index: number of channel to select (index > 0 )     
%
% See also: MTLWCH, MTLWH, MTLRH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995


% if no index available read first channel

if (nargin < 2) index=1; end

iscomplex=0;

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 start,stop,n_directions,comment1, comment2] = mtlrh(filename);

if (file_type == 2) iscomplex=0; end
if (file_type == 3) iscomplex=1; end

[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  channel=-1;
  return
end

pos=info_blocks*256+n_lines*4*(index-1)*(iscomplex+1);

% fseek to channel

status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: mtlrch,  fseek to channel not succesful' );channel=-1; return; end

%fread of channel

ichannel=fread(fid,n_lines*(iscomplex+1),'float32');
if (ichannel<0) error('ERROR: mtlrch, reading of channel not succesful' );channel=-1; end

% creating a complex vector if necessary

if (iscomplex==1) 
   channel(1:n_lines)=ichannel(1:2:n_lines*2)+i*ichannel(2:2:n_lines*2);
 else
   channel=ichannel;
end

% channel=channel';  probably a problem CHK 2/96
channel = channel(:);

status = fclose(fid);
if (status < 0) error('ERROR: mtlrch, close file');channel=-1; return; end

return

function [channel] = MTLreadCH(index,HRTF)
% adapted from MTLRCH"

% MTLRCH reads one channel from a file (type 2 or 3)
%
% [channel]=MTLreadCH(index);
%
% var_string = name of output variable
%
% See also: MTLWCH, MTLWH, MTLreadHDR

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

% if no index available read first channel

if (nargin == 0) index=1; end

file_type = HRTF.filetype;
HRTFfile = [HRTF.directory HRTF.fname];

iscomplex=0;

if (file_type == 2) iscomplex=0; end
if (file_type == 3) iscomplex=1; end

[fid,message] = fopen(HRTFfile,'r','ieee-le');
if fid==-1          % file can't be opened
    disp(message);
    channel=-1;
    return
end

pos=HRTF.info_blocks*256+HRTF.nlines*4*(index-1)*(iscomplex+1);

% fseek to channel
status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: MTLreadCH,  fseek to channel not succesful' );channel=-1; return; end

%fread of channel
ichannel=fread(fid,HRTF.nlines*(iscomplex+1),'float32');
if (ichannel<0) error('ERROR: MTLreadCH, reading of channel not succesful' );channel=-1; end

% creating a complex vector if necessary
if (iscomplex==1)
    channel(1:n_lines)=ichannel(1:2:HRTF.nlines*2)+1i*ichannel(2:2:HRTF.nlines*2);
else
    channel=ichannel;
end

channel = channel(:);

status = fclose(fid);
if (status < 0) error('ERROR: MTLreadCH, close file');channel=-1; return; end
var_string = channel; % sets global variable with name = var_string to value of channel

function [HRTF] = MTLreadHDR(dir,fname)
% Adapted from MTLRH
% MTLreadHDR reads the header of a SIGMAP and MTLSIG file.
%
% filetype: 3 -> frq-domain
%           2 -> time-domain
% info_blocks: number of info blocks
% n_channels: number of channels in file
% n_lines: number of lines per channel
% sampling_rate: sampling rate in kHz
% first_line: value of first line ( in ms for time domain
%                                   in kHz for frq domain )
% last_line: value of last line
% n_directions: number of directions in file
% comment1: string with 80 characters
% comment2: string with 200 characters
%
% See also: MTLWH, MTLRCH, MTLWCH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

%open file for reading
%disp(HRTFfile);
HRTFfile = [dir fname];
[fid,message] = fopen(HRTFfile,'r','ieee-le');
if (fid==-1)         				% file can't be opened
    disp(message);
    filetype=-1;
    return;
end;


fseek(fid,0,'bof');

%test if data are stored in PC-format (EOF=26) or Workstation format(EOF=XX)

fseek(fid,79,'bof');
eof = fread(fid,1,'char');
eof_int = abs(eof);
if (eof_int ~= 26);
    fclose(fid);
    error('chosen file isn''t readable');
    filetype=-1;
    return;
end;

%check filetype

fseek(fid,80,'bof');
filetype = fread(fid,1,'short');
if ( (filetype<2)|(filetype>11) )
    fclose(fid);
    error('choosen file includes unknown filetype');
end;

status=fseek(fid,82,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek info_blocks'); filetype=-1;return, end;
HRTF.info_blocks = fread(fid,1,'short');

status=fseek(fid,84,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek n_channels'); filetype=-1;return, end;
HRTF.nchannels = fread(fid,1,'short');

status=fseek(fid,128,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek n_lines'); filetype=-1;return, end;
HRTF.nlines = fread(fid,1,'short');


%if the filetype is a linear signal (hrir or hrtf)
if (filetype<7)
    status=fseek(fid,132,'bof');
    if (status<0) error('ERROR: MTLreadHDR, fseek sampling_rate'); filetype=-1;return, end;
    HRTF.sampling_rate = fread(fid,1,'float32');
    
    status= fseek(fid,136,'bof');
    if (status<0) error('ERROR: MTLreadHDR, fseek first_line'); filetype=-1;return, end;
    HRTF.first_line = fread(fid,1,'float32');
    
    status=fseek(fid,140,'bof');
    if (status<0) error('ERROR: MTLreadHDR, fseek last_lines'); filetype=-1;return, end;
    HRTF.last_line = fread(fid,1,'float32');
end;


% if the file is sigmap format it does not contain 'n_directions' field
% assume that it has one direction

status=fseek(fid,220,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek n_directions'); filetype=-1;return, end;
HRTF.ndirections = fread(fid,1,'short');
if (HRTF.ndirections==0) HRTF.ndirections=1; end;

status=fseek(fid,0,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek comment1'); filetype=-1;return, end;
comment1 = setstr(fread(fid,80,'uchar')) ;
HRTF.comment1=comment1';

status=fseek(fid,222,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek comment2'); filetype=-1;return, end;
comment2 = setstr(fread(fid,200,'uchar'));
HRTF.comment2=comment2';

HRTF.filetype = filetype;
fclose(fid);
return;
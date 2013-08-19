function  mtlwh(filename,filetype,info_blocks,...
          n_channels,n_lines, sampling_rate,...
          first_line,last_line,n_directions,...
          comment1, comment2)

%function  mtlwh(filename,filetype,info_blocks,...
%          n_channels,n_lines, sampling_rate,...
%          first_line,last_line,n_directions,...
%          comment1, comment2)
%
% MTLWH writes the header of a SIGMAP and MTLSIG file.
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
% See also: MTLRH, MTLRCH, MTLWCH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

%open file for reading and appending new data

% [fid,message] = fopen(filename,'a+','ieee-le');
if exist(filename)==2		% file exists
   [fid,message] = fopen(filename,'r+','l');
else				% file does not exist
   [fid,message] = fopen(filename,'w+','l');
end

if fid==-1          % file can't be opened
  disp(message);
  return;
end

%erase old info blocks or create new info blocks
%commented out Dec 4, 1995
 status=fseek(fid,0,'bof');
 if (status < 0) error('ERROR: mtlwh, erase old info blocks, fseek'); return; end
 d=zeros(1,256*info_blocks);
 count = fwrite(fid,d,'uchar');
 if (count < (256*2)) error('ERROR: mtlwh, erase old info blocks, fwrite'); return; end
%end of Dec 4 mod

% write comment
if (filetype<=3)
   status=fseek(fid,0,'bof');
   if (status < 0) error('ERROR: mtlwh, write comment, fseek'); return; end
   n = length (comment1);
   comment1=[comment1 zeros(1,79-n)];
   count=fwrite(fid,comment1,'uchar');
   if (count < length(comment1)) error('ERROR: mtlwh, write comment, fwrite'); return; end
end

%write EOF
status=fseek(fid,79,'bof');
if (status < 0) error('ERROR: mtlwh, write EOF, fseek'); return; end
eof=26;
count = fwrite(fid,eof,'char');
if (count < length(eof)) error('ERROR: mtlwh, write EOF, fwrite'); return; end


%write filetype
status=fseek(fid,80,'bof');
if (status < 0) error('ERROR: mtlwh, write filetype, fseek'); return; end
count = fwrite(fid,filetype,'short');
if (count < length(filetype)) error('ERROR: mtlwh, write filetype, fwrite'); return; end


%write info_blocks
status=fseek(fid,82,'bof');
if (status < 0) error('ERROR: mtlwh, write info_blocks, fseek'); return; end
count = fwrite(fid,info_blocks,'short');
if (count < length(info_blocks)) error('ERROR: mtlwh, write info_blocks, fwrite'); return; end


%write n_channels
status=fseek(fid,84,'bof');
if (status < 0) error('ERROR: mtlwh, write n_channels, fseek'); return; end
count = fwrite(fid,n_channels,'short');
if (count < length(n_channels)) error('ERROR: mtlwh, write n_channels, fwrite'); return; end


%write n_lines
status=fseek(fid,128,'bof');
if (status < 0) error('ERROR: mtlwh, write n_lines, fseek'); return; end
count = fwrite(fid,n_lines,'long');
if (count < length(n_lines)) error('ERROR: mtlwh, write n_lines, fwrite'); return; end


if (filetype<7)
   %write sampling rate 
  status=fseek(fid,132,'bof');
  if (status < 0) error('ERROR: mtlwh, write sampling_rate, fseek'); return; end
  count = fwrite(fid,sampling_rate,'float32');
  if (count < length(sampling_rate)) error('ERROR: mtlwh, write sampling_rate, fwrite'); return; end


  %write first_line
  status=fseek(fid,136,'bof');
  if (status < 0) error('ERROR: mtlwh, write first_line, fseek'); return; end
  count = fwrite(fid,first_line,'float32');
  if (count < length(first_line)) error('ERROR: mtlwh, write first_line, fwrite'); return; end


  %write last_line

  status=fseek(fid,140,'bof');
  if (status < 0) error('ERROR: mtlwh, write last_line, fseek'); return; end
  count = fwrite(fid,last_line,'float32');
  if (count < length(last_line)) error('ERROR: mtlwh, write last_line, fwrite'); return; end
end


%write n_directions
status=fseek(fid,220,'bof');
if (status < 0) error('ERROR: mtlwh, write n_directions, fseek'); return; end
count = fwrite(fid,n_directions,'short');
if (count < length(n_directions)) error('ERROR: mtlwh, write n_directions, fwrite'); return; end

%write comment2
status=fseek(fid,222,'bof');
if (status < 0) error('ERROR: mtlwh, write comment2, fseek'); return; end
comment1=[comment1 zeros(1,200-n)];
count= fwrite(fid,comment2,'uchar');
if (count < length(comment2)) error('ERROR: mtlwh, write comment2, fwrite'); return; end

%write something to fill end of info block
status= fseek(fid,511,'bof');
if (status < 0) error('ERROR: mtlwh, write something, fseek'); return; end
fwrite(fid,'a','char');
if (count < 1) error('ERROR: mtlwh, write something, fwrite'); return; end


%close file
status=fclose(fid);
if (status < 0) error('ERROR: mtlwh, close file'); return; end
function status = mtlwch(filename,channel, index)

% status = mtlwch(filename,channel, index)
%
% MTLWCH writes one channel to a file (type 2 or 3)
%
%         filename: string with filename
%         channel: vector with data
%         index: number of channel (index > 0)       
%         
% See also: MTLRCH, MTLRH, MTLWH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

iscomplex=0;

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 start,stop,n_directions,comment1, comment2] = mtlrh(filename);


if (file_type == 2) iscomplex=0; end
if (file_type == 3) iscomplex=1; end

%[fid,msg]=fopen(filename,'a+','ieee-le');
if exist(filename)==2		% file exists
   [fid,message] = fopen(filename,'r+','l');
else				% file does not exist
   [fid,message] = fopen(filename,'w+','l');
end;


if fid==-1          % file can't be opened
  disp(message);
  status=-1;
  return;
end

status=fseek(fid,0,'eof');
if (status<0) error('ERROR: mtlwch,  fseek to end not succesful' ); status=-1;return; end
end_pos=ftell(fid);
if (end_pos < 0) error('ERROR: mtlwch,  ftell at end not succesful' ); status=-1;return; end

pos=info_blocks*256+n_lines*4*(index-1)*(iscomplex+1);
if (pos > end_pos) error('ERROR: pos > end_pos'); status=-1; return; end

status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: mtlwch,  fseek to channel position not succesful' ); status=-1;return;
end

if(iscomplex==1)
 ochannel(1,:)=(real(channel))';					%CHK
 ochannel(2,:)=(imag(channel))';					%CHK
 %ochannel(1,:)=real(channel);
 %ochannel(2,:)=imag(channel);
 ochannel=reshape(ochannel,1,n_lines*2);
 %disp(['data complex']);
else
 %disp(['data not complex']);
 ochannel=channel';
end

count=fwrite(fid,ochannel','float32');
disp(['wrote ' num2str(count) ' points']);
end_pos=ftell(fid);

if (count<(n_lines*(iscomplex+1))) 
% if (ichannel<0) error('ERROR: mtlwch, writing of channel not successful' );status=-1; return; end;
 error('ERROR: mtlwch, writing of channel not successful' );   % CHK 2/96
 status=fclose(fid);
 if (status < 0) error('ERROR: mtlwch, close file');status=-1; return; end
 return;
end
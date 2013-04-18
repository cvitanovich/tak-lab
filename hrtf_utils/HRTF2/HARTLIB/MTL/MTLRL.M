function output=mtlrl(filename, line);

% MTLLR reads one line across from a file (type 2 or 3)
%
% [channel]=mtlrl(filename, line);
%
%         filename: string with filename
%         line: number of channel to select (line > 0 )     
%
% See also: MTLWCH, MTLWH, MTLRH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 20.4.1995




% if no index available read first line

if (nargin < 2) line=0; end;
iscomplex=0;

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 start,stop,n_directions,comment1, comment2] = mtlrh(filename);

if (line > n_lines) error('ERROR: mtlrh, line > n_lines');output=-1; return; end;

if (file_type == 2) iscomplex=0; end;
if (file_type == 3) iscomplex=1; end;


[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  channel=-1;
  return;
end;


output=zeros(1,n_directions);


for a_dir=1:n_directions
	 pos=info_blocks*256+(n_lines*(a_dir-1)+line)*(iscomplex+1)*4;
         
	 % fseek to line

	 status=fseek(fid,pos,'bof');
	 if (status<0) error('ERROR: mtlrh,  fseek to channel not succesful' );channel=-1; return; end;


	 %fread of channel

	 ichannel=fread(fid,iscomplex+1,'float32');
%if (ichannel<0) error('ERROR: mtlrh, reading of channel not succesful' );channel=-1; end;


	% creating a complex vector if necessary

	if (iscomplex==1) 
	   channel=ichannel(1)+i*ichannel(2);
 	else
	   channel=ichannel;
	end

	output(a_dir)=channel;

end;

status = fclose(fid);
if (status < 0) error('ERROR: mtlrh, close file');output=-1; return; end;

return;

end;

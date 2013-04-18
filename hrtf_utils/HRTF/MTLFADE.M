function mtlfade(infile,outfile,index)

% function mtlfade(infile,outfile,indices)
% notice parameters (index) for where to fade
% these have been changed 1/22/98
% and will change again for the monkey

if (nargin < 3)
%   index = [171, 271, 839, 936];     % old owl params
   index = [103, 171, 839, 936];     % new owl params 1/22/98 = 1496, 2500, 12287, 13710 Hz
%   index = [1, 30, 1475, 1639];     % monkey params 1/22/98 = 0, 330, 18002, 20005 Hz
end

if (nargin<2) 
disp(['outfile not specified']);
  return
end;

[filetype,info_blocks,...
 n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,...
 comment1, comment2] = mtlrh(infile);

if (filetype ~= 3) 
  disp(['File: ' infile ' has wrong type']);
  return;
end;


binaural=(n_channels/n_directions)-1;
index=0;

%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile,n_directions);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
      delay=mtlrdel2(infile);
   end;
end;



n_dir=n_directions;
if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256)
  info_blocks=2+n_blocks;
end;


filetype=3;


mtlwh(outfile,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);



%[fid,message] = fopen(outfile,'a','ieee-le');
if exist(outfile)==2		% file exists
   [fid,message] = fopen(outfile,'r+','l');
else				% file does not exist
   [fid,message] = fopen(outfile,'w+','l');
end;


if fid==-1          % file can't be opened
  disp(message);
  return;
end;

 

%if direction_matrix
if ((n_blocks) > 0)
 status=mtlwdir(outfile,direction_matrix);
 if(length(delay) > 0) mtlwdel2(outfile, delay); end;
end


epsilon=1E-8;

for index=1:n_channels
      clc
      disp(['MTLFADE: Calculating channel ', num2str(index), ' of ', num2str(n_channels)] );
      channel = mtlrch(infile,index);
      channel=channel(:);
      %channel=fade(channel,171,271,839,936);   % used to be 939 CHK 3/18/96
      channel = fade(channel,index(1), index(2), index(3), index(4));   % changed 1/22/98
      p=find(channel==0); channel(p)=ones(size(p))*epsilon;
      status = mtlwch(outfile,channel,index);
      end;
return;
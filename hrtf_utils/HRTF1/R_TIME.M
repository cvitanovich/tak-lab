function [file_format,first_channel,second_channel,direction_matrix,time_delay,...
          n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,...
          end_of_signal,xaxis,comment] = r_time(fid,filetype)

% function [file_format,first_channel,second_channel,direction_matrix,time_delay,...
%           n_directions,n_channels,n_lines,sampling_rate,...
%           begin_of_signal,end_of_signal,xaxis] = r_time(fid,filetype)
%
% R_TIME Reads the specific time signal attributes.
%
%        meaning of the input arguments:
%         - fid: internal file identification number of the with the function
%                'read' opened file
%         - filetype: number of identification for the signal type
%
%        meaning of the output arguments:
%         - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                        labeling
%         - first_channel: matrix with the absolute magnitude of the signals
%                          first channel. Each row represents one direction,
%                          each column one sampling line.
%         - second_channel: matrix with the absolute magnitude of the signals
%                           second channel. Each row represents one 
%                           direction, each column one sampling line. If
%                           the file includes only one channel, this matrix
%                           is empty.
%         - direction_matrix: matrix including the elevation and azimuth
%                             angle of the different directions. The first
%                             row includes the elevation, the second row
%                             the azimuth. Each column represents a
%                             direction.
%         - time_delay: matrix including the initial time delays of the
%                       impulse responses. Each column represents one direc-
%                       tion, each row one channel.
%         - n_directions: number of directions
%         - n_channels: number of channels
%         - n_lines: number of lines
%         - sampling_rate
%         - begin_of_signal
%         - end_of_signal
%         - xaxis: vector including the points of sampling
%         - comment: commentary included in MTLSIG and SIGMAP files.
%
%        For further information see: MTLSIG manuals, SIGMAP manuals, READ
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


fseek(fid,82,'bof');
info_blocks = fread(fid,1,'short');
file_channels = fread(fid,1,'short');
fseek(fid,128,'bof');
n_lines = fread(fid,1,'int32');
if (filetype<7)
  sampling_rate = fread(fid,1,'float32');
else
  fseek(fid,136,'bof');
end;
begin_of_signal = fread(fid,1,'float32');
end_of_signal = fread(fid,1,'float32');

fseek(fid,220,'bof');
n_directions = fread(fid,1,'short');
if n_directions==0
 n_directions = file_channels;
 file_format = 'SIGMAP';
 disp('SIGMAP format loaded');
 begin_of_signal = 0;
 end_of_signal = 1 / sampling_rate * (n_lines-1);
else
 file_format = 'MTLSIG';
 disp('MTLSIG format loaded');
end;

if file_format=='MTLSIG'
  comment = setstr(fread(fid,200,'uchar'));
else
  fseek(fid,0,'bof');
  comment = setstr(fread(fid,80,'uchar'));    
end;

n_channels = file_channels / n_directions;

disp(['file type: ' num2str(filetype)]);
disp(['number of lines: ' num2str(n_lines)]);
disp(['number of channels: ' num2str(n_channels)]);
disp(['number of directions: ' num2str(n_directions)]);
disp('comment: ');disp(comment');
if (filetype<7)
  disp(['sampling rate: ',num2str(sampling_rate),' Hz']);
end;
disp(['begin: ',num2str(begin_of_signal),' sec']);
disp(['end: ',num2str(end_of_signal),' sec']);

if file_format=='MTLSIG'
  fseek(fid,512,'bof');
  direction_matrix = fread(fid,2*n_directions,'float32');
  direction_matrix = reshape(direction_matrix ,2,n_directions);
  time_delay = fread(fid,n_channels*n_directions,'float32');
  time_delay = reshape(time_delay,n_channels,n_directions);
end;

xaxis = linspace(begin_of_signal,end_of_signal,n_lines);

fseek(fid,info_blocks*256,'bof');

for i=1:n_directions
   if n_channels==1
     row = fread(fid,n_lines,'float32');
     row = row';
     first_channel = [first_channel;row];
   else
     row = fread(fid,n_lines,'float32');
     row = row';
     first_channel = [first_channel;row];
     row = fread(fid,n_lines,'float32');
     row = row';
     second_channel = [second_channel;row];
   end
end

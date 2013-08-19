function [file_format,first_channel,second_channel,first_phase,second_phase,...
          direction_matrix,n_directions,n_channels,n_lines,sampling_rate,...
          begin_of_signal,end_of_signal,xaxis,comment] = r_frq(fid,filetype)

% function [file_format,first_channel,second_channel,first_phase,second_phase,direction_matrix,...
%           n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,...
%           end_of_signal,xaxis,comment] = r_frq(fid,filetype)
%
% R_FRQ Reads the specific frequency signal attributes. 
%
%       meaning of the input arguments:
%        - fid: internal file identification number of the with the function
%               'read' opened file
%        - filetype: number of identification for the signal type
%
%       meaning of the output arguments:
%        - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                       labeling
%        - first_channel: matrix with the absolute magnitude of the signals
%                         first channel. Each row represents one direction,
%                         each column one sampling line. Reading interaural
%                         differences, 'first_channel' includes the interaural
%                         level differences.
%        - second_channel: matrix with the absolute magnitude of the signals
%                          second channel. Each row represents one direction,
%                          each column one sampling line. If the file includes
%                          only one channel, this matrix is empty.
%        - first_phase: phase, phase delay or group delay of the first
%                       channel.  Each row represents one direction, each
%                       column one sampling line. Reading interaural
%                       differences, 'first_phase' includes the interaural
%                       time differences.
%        - second_phase: phase, phase delay or group delay of the second
%                        channel, if the file includes two channels. Else
%                        the variable is empty. Each row represents one
%                        direction, each column one sampling line.
%        - direction_matrix: matrix including the elevation and azimuth
%                            angle of the different directions. The first
%                            row includes the elevation, the second row
%                            the azimuth. Each column represents a
%                            direction.
%        - n_directions: number of directions
%        - n_channels: number of channels
%        - n_lines: number of lines
%        - sampling_rate
%        - begin_of_signal
%        - end_of_signal
%        - xaxis: vector including the points of sampling
%        - comment: commentary included in MTLSIG and SIGMAP files.
%
%       For further information see: MTLSIG manuals, SIGMAP manuals, READ 
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


info_blocks = fread(fid,1,'short');
file_channels = fread(fid,1,'short');
fseek(fid,128,'bof');
n_lines = fread(fid,1,'long');
if (filetype<7)|(filetype==11)
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
  n_lines = n_lines / 2 + 1;
     % n_lines sampled but only half of them stored because of periodicity
     % and real part at Nyquist-frequency
  file_format = 'SIGMAP';
  disp('SIGMAP format loaded');
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
    % number_of_channels==1, when reading a SIGMAP file
    % each channel taken as a different direction

disp(['file type: ' num2str(filetype)]);
disp(['number of lines: ' num2str(n_lines)]);
disp(['number of channels: ' num2str(n_channels)]);
disp(['number of directions: ' num2str(n_directions)]);
disp('comment: ');disp(comment');
if (filetype<7)|(filetype==11)
  disp(['sampling rate: ',num2str(sampling_rate),' Hz']);
end;
if ( (filetype>=3)&(filetype<=6) )
  disp(['begin: ',num2str(begin_of_signal),' Hz']);
  disp(['end: ',num2str(end_of_signal),' Hz']);
elseif filetype==7
  disp(['begin: ',num2str(begin_of_signal),' mel']);
  disp(['end: ',num2str(end_of_signal),' mel']);
elseif filetype==8
  disp(['begin: ',num2str(begin_of_signal),' bark']);
  disp(['end: ',num2str(end_of_signal),' bark']);
elseif filetype==9
  disp(['begin: ',num2str(begin_of_signal),' ERB']);
  disp(['end: ',num2str(end_of_signal),' ERB']);
elseif filetype==10
  disp(['begin: ',num2str(begin_of_signal),' log10( f / Hz )']);
  disp(['end: ',num2str(end_of_signal),' log10( f / Hz)']);
elseif filetype==11
  disp(['begin: ',num2str(begin_of_signal),' Hz']);
  disp(['end: ',num2str(end_of_signal),' Hz']);
end;

if file_format=='MTLSIG'
  fseek(fid,512,'bof');
  direction_matrix = fread(fid,2*n_directions,'float32');
  direction_matrix = reshape(direction_matrix,2,n_directions);
end;

if (filetype>=7)&(filetype<=11)
  xaxis = fread(fid,n_lines,'float32');
  xaxis = xaxis';
else
  xaxis = linspace(begin_of_signal,end_of_signal,n_lines);
end;

fseek(fid,info_blocks*256,'bof');

if ( (filetype==3)|(filetype==7)|(filetype==8)|(filetype==9)|(filetype==10) )
          % stored with real- and imaginary part
  if n_channels==1
    for r=1:n_directions
      disp(['R_FRQ: Reading direction # ' num2str(r)]);
      if file_format=='MTLSIG'
        one_direct = fread(fid,n_lines*2,'float32');
        one_direct = reshape(one_direct,2,n_lines);
        r_part = one_direct(1,:);
        i_part = one_direct(2,:);
        channel = r_part+j*i_part;
      else                    % SIGMAP format
        one_direct = fread(fid,(n_lines-1)*2,'float32');
        one_direct = reshape(one_direct,2,(n_lines-1));
        r_part = one_direct(1,:);
        i_part = one_direct(2,:);
        channel = r_part+j*i_part;
        channel(n_lines) = i_part(1)+j*0;  % Nyquist part
        channel(1) = real(channel(1));
      end;
      first_channel = [first_channel;channel];
    end %for
  else                                     % number_of_channels==2
    for r=1:n_directions
      disp(['R_FRQ: Reading direction # ' num2str(r)]);
      one_direct = fread(fid,n_lines*2,'float32');     % reading first channel
      one_direct = reshape(one_direct,2,n_lines);
      r_part = one_direct(1,:);
      i_part = one_direct(2,:);
      channel = r_part+j*i_part;
      first_channel = [first_channel;channel];

      one_direct = fread(fid,n_lines*2,'float32');     % reading second channel
      one_direct = reshape(one_direct,2,n_lines);
      r_part = one_direct(1,:);
      i_part = one_direct(2,:);
      channel = r_part+j*i_part;
      second_channel = [second_channel;channel];
    end % for (n_directions)
  end % if
 % end of filetype 3,7,8,9,10
else
 % stored with magnitude and phase (filetype 3) or
 % interaural level difference and time difference (filetype 11)
 if n_channels==1
    for r=1:n_directions
      disp(['R_FRQ: Reading direction # ' num2str(r)]);
      one_direct = fread(fid,n_lines*2,'float32');
      one_direct = reshape(one_direct,2,n_lines);
      magnitude = one_direct(1,:);
      phase = one_direct(2,:);
      first_channel = [first_channel;magnitude];
      first_phase = [first_phase;phase];
    end %for
 else                                     % number_of_channels==2
    for r=1:n_directions
      disp(['R_FRQ: Reading direction # ' num2str(r)]);
      one_direct = fread(fid,n_lines*2,'float32');     % reading first channel
      one_direct = reshape(one_direct,2,n_lines);
      magnitude = one_direct(1,:);
      phase = one_direct(2,:);
      first_channel = [first_channel;magnitude];
      first_phase = [first_phase;phase];

      one_direct = fread(fid,n_lines*2,'float32');     % reading second channel
      one_direct = reshape(one_direct,2,n_lines);
      magnitude = one_direct(1,:);
      phase = one_direct(2,:);
      second_channel = [second_channel;magnitude];
      second_phase = [second_phase;phase];
    end % for (n_directions)
  end % if
end % end of filetype 4,5,6


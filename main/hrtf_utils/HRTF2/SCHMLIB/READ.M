function [file_format,filetype,first_channel,second_channel,first_phase,...
          second_phase,direction_matrix,time_delay,n_directions,n_channels,...
          n_lines,sampling_rate,begin_of_signal,end_of_signal,...
          xaxis,comment] = read(filename)

% function [file_format,filetype,first_channel,second_channel,first_phase,second_phase,...
%           direction_matrix,time_delay,n_directions,n_channels,n_lines,...
%           sampling_rate,begin_of_signal,end_of_signal,xaxis,comment]...
%           = read(filename)
%
% READ Reads the catalog of hrtf or hrir signals,if the file format is the
%      MTLSIG - format. If the fileformat is the SIGMAP - format the impulse
%      responses or transfer functions are read, with out any spatial order.
%      For a SIGMAP - file all signals are treated as signals of the same channel
%      with different directions. For different layouts using the provided
%      display functions, the output argument 'file_format' distinguishes
%      between the two file formats. Attributes that are not provided by
%      SIGMAP files return as empty elements. 'n_directions' includes in this
%      case the number of different signals and 'n_channels' is one.
%
%      The input argument 'filename' includes the name of the MTLSIG file,
%      that should be opened.
%
%      The following variables are returned:
%       - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                      labeling
%       - filetype: distinguishes between impulse responses, the different
%                   format types of transfer functions and interaural differences
%       - first_channel: matrix with the absolute magnitude of the signals
%                        first channel. Each row represents one direction,
%                        each column one sampling line. Reading interaural
%                        differences, 'first_channel' includes the interaural
%                        level differences.
%       - second_channel: matrix with the absolute magnitude of the signals
%                         second channel. Each row represents one direction,
%                         each column one sampling line. If the file includes
%                         only one channel, this matrix is empty.
%       - first_phase: phase, phase delay or group delay of the first
%                      channel.  Each row represents one direction, each
%                      column one sampling line. Reading interaural
%                      differences, 'first_phase' includes the interaural
%                      time differences.
%       - second_phase: phase, phase delay or group delay of the second
%                       channel, if the file includes two channels. Else
%                       the variable is empty. Each row represents one
%                       direction, each column one sampling line.
%       - direction_matrix: matrix including the elevation and azimuth
%                           angle of the different directions. The first
%                           row includes the elevation, the second row
%                           the azimuth. Each column represents a
%                           direction.
%       - time_delay: matrix including the initial time delays of the
%                     impulse responses. Each column represents one direc-
%                     tion, each row one channel.
%       - n_directions: number of directions
%       - n_channels: number of channels
%       - n_lines: number of lines for each direction and each channel
%       - sampling_rate
%       - begin_of_signal
%       - end_of_signal
%       - xaxis: vector including the points of sampling
%       - comment: commentary included in MTLSIG and SIGMAP files.
%
%      Important !
%      This function may return empty elements, if the signal type or file format
%      does not provide such a certain attribute. For example: time signals
%      have no phase and transfer functions have not any time_delay. Which
%      elements a signal includes is written in the SIGNAL MANUAL PAGES.
%
%      For further information see: MTLSIG manuals, SIGMAP manuals
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


disp(['file: ',filename]);
[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

fseek(fid,79,'bof');
eof = fread(fid,1,'char');
eof_int = abs(eof);
if (eof_int ~= 26)
  fclose(fid);
  error('choosen file isn''t readable');
end;

fseek(fid,80,'bof');
filetype = fread(fid,1,'short');
if ( (filetype<2)|(filetype>11) )
  fclose(fid);
  error('choosen file includes unknown filetype');
end;
 
if (filetype==2)                       % reading a time signal
  [file_format,first_channel,second_channel,direction_matrix,time_delay,...
   n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,...
   end_of_signal,xaxis,comment] = r_time(fid,filetype);
else                                    % reading a frq signal
  [file_format,first_channel,second_channel,first_phase,second_phase,...
   direction_matrix,n_directions,n_channels,n_lines,sampling_rate,...
   begin_of_signal,end_of_signal,xaxis,comment] = r_frq(fid,filetype);
end

fclose(fid);

function meshplot(file_format,filetype,filename,first_channel,second_channel,...
                  first_phase,second_phase,direction_matrix,...
                  n_directions,n_channels,n_lines,xaxis,scaling)

% function meshplot(file_format,filetype,filename,first_channel,second_channel,...
%                   first_phase,second_phase,direction_matrix,...
%                   n_directions,n_channels,n_lines,xaxis,scaling)
%
% MESHPLOT generates a meshplot of the magnitude for each channel. For
%          frequency signals the phase for each channel is plotted, too.
%
%          meaning of the input arguments:
%          - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                         labeling
%          - filetype: number of identification for the signal type
%          - filename: name of the file
%          - first_channel: matrix with the absolute magnitude of the
%                           signals first channel
%          - second_channel: matrix with the absolute magnitude of the signals
%                            second channel. If the file includes only one
%                            channel, this matrix is empty.
%          - first_phase: matrix with the phase of the first channel. If the
%                         file includes a time signal, the matrix is empty.
%          - second_phase: matrix with the phase of the second channel. If
%                          the file includes a time signal or a frequency
%                          signal with only one channel, this matrix is
%                          empty.
%          - direction_matrix: matrix including the elevation and azimuth angle
%                              of the different directions. The first row
%                              includes the elevation, the second row the
%                              azimuth. Each column represents a direction.
%          - n_directions: number of directions
%          - n_channels: number of channels
%          - n_lines: number of lines
%          - xaxis: vector including the points of sampling
%          - scaling: must be 'lin' or 'log' and defines the scaling of the
%                     frequency axis
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

if filetype==2    % time signal
  if n_channels==1
    signal = first_channel;
    meshing;
    title(['impulse response  ',filename]);
  else
    signal = first_channel;
    meshing;
    title(['impulse resonse first channel  ',filename]);
    signal =second_channel;
    meshing;
    title(['impulse response second channel  ',filename]);
  end
else                          % frequency signal
  if n_channels==1
    signal = abs2db(first_channel);
    z_status = 0;  % status for zlabel
    meshing;
    title(['transfer function: magnitude first channel  ',filename]);
    if ( (filetype>=4)&(filetype<=6) )
      signal = first_phase;
    else
      signal = [];
      for i=1:n_directions
       sig = 180/pi*unwrap(angle(first_channel(i,:)));
       signal = [signal;sig];
      end;
    end;
    z_status = 1;  % status for zlabel
    meshing;
    title(['transfer function: phase first channel  ',filename]);
  else
 % plotting first channel

    signal = abs2db(first_channel);
    z_status = 0;  % status for zlabel
    meshing;
    title(['transfer function: magnitude first channel  ',filename]);
    if ( (filetype>=4)&(filetype<=6) )
      signal = first_phase;
    else
      signal = [];
      for i=1:n_directions
       sig = 180/pi*unwrap(angle(first_channel(i,:)));
       signal = [signal;sig];
      end;
    end;
    z_status = 1;  % status for zlabel
    meshing;
    title(['transfer function: phase first channel  ',filename]);

 % plotting second channel

    signal = abs2db(second_channel);
    z_status = 0;
    meshing;
    title(['transfer function: magnitude second channel  ',filename]);
    if ( (filetype>=4)&(filetype<=6) )
      signal = second_phase;
    else
      signal = [];
      for i=1:n_directions
       sig = 180/pi*unwrap(angle(second_channel(i,:)));
       signal = [signal;sig];
      end;
    end;
    z_status = 1;
    meshing;
    title(['transfer function: phase second channel  ',filename]); 
  end % end n_channels
end  % end time- or frequency signal


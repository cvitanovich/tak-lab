function contplot(file_format,filetype,filename,first_channel,second_channel,...
                  first_phase,second_phase,direction_matrix,n_directions,...
                  n_channels,n_lines,xaxis,scaling)

% function contplot(file_format,filetype,filename,first_channel,second_channel,...
%                   first_phase,second_phase,direction_matrix,n_directions,...
%                   n_channels,n_lines,xaxis,scaling)
%
% CONTPLOT plots for frequency signals a contour plot of each channels
%          magnitude and phase. The datas are resampled for smoothing the
%          contour lines. If a contour plot of the magnitude is done, the 
%          contour lines are selected by choosing the level difference.
%          Contour plots of the phase are done with 10 different level lines.
%          An axis-control allows the zooming of a contour plot. The contour
%          lines are labeled manual. Only frequency signals can be plotted
%          with contour lines.
%
%          meaning of the input arguments:
%           - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                          labeling
%           - filetype: number of identification for the signal type
%           - filename: name of the file
%           - first_channel: matrix with the absolute magnitude of the signals
%                            first channel
%           - second_channel: matrix with the absolute magnitude of the
%                             signals second channel. If the file includes
%                             only one channel, this matrix is empty.
%           - first_phase: phase, phase delay or group delay of the first
%                          channel
%           - second_phase: phase, phase delay or group delay of the second
%                           channel, if the file includes two channels. Else
%                           the variable is empty.
%           - direction_matrix: matrix including the elevation and azimuth
%                               angle of the different directions. The first
%                               row includes the elevation, the second row
%                               the azimuth. Each column represents a
%                               direction.
%           - n_directions: number of directions
%           - n_channels: number of channels
%           - n_lines: number of lines
%           - xaxis: vector including the points of sampling
%           - scaling: must be 'lin' or 'log' and defines the scaling of the
%                      frequency axis
%
%          For further informations see: MTLSIG - manuals
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

if filetype==2    			% time signal
  disp('contour plot for time signals not availabel');
else                         		% frequency signal
  if n_channels==1
    signal = abs2db(first_channel);
    status = 0; % status for zlabel
    cont_ing;
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
    status = 1; 			% status for zlabel
    cont_ing;
    title(['transfer function: phase first channel  ',filename]);
  else					% plotting first channel
    signal = abs2db(first_channel);
    status = 0;				% status for zlabel
    cont_ing;
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
    status = 1; 			% status for zlabel
    cont_ing;
    title(['transfer function: phase first channel  ',filename]);

 					% plotting second channel

    signal = abs2db(second_channel);
    status = 0; 			% status for zlabel
    cont_ing;
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
    status = 1; 			% status for zlabel
    cont_ing;
    title(['transfer function: phase second channel  ',filename]); 
  end 					% end n_channels
end  					% end time- or frequency signal


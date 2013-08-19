function t_layout(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,...
                  direction_matrix,xaxis,n_lines)

% function t_layout(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,...
%                   direction_matrix,xaxis,n_lines)
%
% T_LAYOUT generates a layout for plots of impulse responses. Different
%          layouts are available, depending on the different file formats.
%
%          The following attributes are provided:
%
%          - controls the correct range for the x-axis
%          - labels the axes
%          - displayes title
%          - displayes the directions
%
%          meaning of the input arguments:
%          - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                         labeling
%          - filetype: number of identification for the signal type
%          - filename: name of the file
%          - n_o_subplot: number of subplot in the actual figure (only the
%                         first subplot has a title)
%          - headline: title of the first subplot
%          - n_o_direction: number of direction; index about all plotted
%                           signals
%          - direction_matrix: matrix including the elevation and azimuth
%                              angle of the different directions. The first
%                              row includes the elevation, the second row
%                              the azimuth. Each column represents a
%                              direction.
%          - xaxis: vector including the points of sampling
%          - n_lines: number of lines
%
%          For further information see: magn_plot
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


% axis changing
axissetting = axis;
axissetting(1) = xaxis(1);
axissetting(2) = xaxis(n_lines);
axis(axissetting);


% title and labels
if filetype==2                            % time signal
  if n_o_subplot==1
    title([headline,'   ',filename]);
  end
  ylabel('h(t)');
  xlabel('t / sec');
else                                       % frequency signal
 error('choose another layout function');
end

% direction display
if file_format=='MTLSIG'
  dirlabel(filetype,n_o_direction,direction_matrix);
end;
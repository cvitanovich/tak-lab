function delayout(filetype,filename,n_o_subplot,headline,n_o_direction,...
                  direction_matrix,xaxis,n_lines)

% function delayout(filetype,filename,n_o_subplot,headline,n_o_direction,...
%                   direction_matrix,xaxis,n_lines)
%
% DELAYOUT generates a layout for 'delay plots'
%          - controls the correct range for the x-axis
%          - labels the axes
%          - displayes title
%          - displayes the directions
%
%          meaning of the input arguments:
%          - filetype: number of identification for the signal type
%          - filename: name of the file
%          - n_o_subplot: number of subplot in the actual figure (only the
%                         first subplot has a title)
%          - headline: title of the first subplot
%          - n_o_direction: number of direction; is an index about all plotted
%                           signals
%          - direction_matrix: matrix including the elevation and azimuth angle
%                              of the different directions. The first row
%                              includes the elevation, the second row the
%                              azimuth. Each column represents a direction.
%          - xaxis: vector including the points of sampling
%          - n_lines: number of lines
%
%          For further information see: delay
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
if ( (filetype==5)|(filetype==6)|(filetype==4) )          % bode signal
 if n_o_subplot==1
   title([headline,filename]);
 end;
 if filetype==5 
   ylabel('B(f) / f');
 else ylabel('B''(f) / f''');

 xlabel('f / Hz');

 end;
else
 error('choose another layout function');
end

% direction display
dirlabel(filetype,n_o_direction,direction_matrix);

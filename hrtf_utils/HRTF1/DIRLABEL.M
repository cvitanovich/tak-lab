function dirlabel(filetype,n_o_direction,direction_matrix)

% function dirlabel(filetype,n_o_direction,direction_matrix)
%
% DIRLABEL labels the plots with the direction belonging to it.
%          For impulse responses the label is placed in the upper right
%          corner, for transfer functions in the lower left corner.
%
%          meaning of the input arguments:
%          - filetype: number of identification for the signal type
%          - n_o_direction: number of direction; is an index about all signals
%          - direction_matrix: matrix including the elevation and azimuth angle
%                              of the different directions. The first row
%                              includes the elevation, the second row the
%                              azimuth. Each column represents a direction.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

elev = num2str(direction_matrix(1,n_o_direction));
azi = num2str(direction_matrix(2,n_o_direction));
e1 = text('String','elev: ');
e2 = text('String',elev);
a1 = text('String','azi: ');
a2 = text('String',azi);
if filetype==2
  set(e1,'Units','normalized','Position', [.85 .85]);
  set(e2,'Units','normalized','Position', [.93 .85]);
  set(a1,'Units','normalized','Position', [.85 .6]);
  set(a2,'Units','normalized','Position', [.93 .6]);
else
  set(e1,'Units','normalized','Position', [.02 .4]);
  set(e2,'Units','normalized','Position', [.1 .4]);
  set(a1,'Units','normalized','Position', [.02 .15]);
  set(a2,'Units','normalized','Position', [.1 .15]);
end

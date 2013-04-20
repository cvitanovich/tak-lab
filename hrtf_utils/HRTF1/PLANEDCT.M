function single_plane = planedct(direction_matrix)

% function single_plane = planedct(direction_matrix)
%
% PLANEDCT detects if the direction matrix includes a single horizontal
%          plane. Then 'single_plane' includes the string 'elev', else the
%          string 'none' is returned.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


min_elev = min(direction_matrix(1,:));
max_elev = max(direction_matrix(1,:));

if min_elev==max_elev
 single_plane = 'elev';
else
 single_plane = 'none';
end;
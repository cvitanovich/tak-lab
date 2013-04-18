function [ytic,lab,ang1] = y_axis(direction_matrix,n_directions)

% function [ytic,lab,ang1] = y_axis(direction_matrix,n_directions)
%
% Y_AXIS initializes a vector 'ytic' with the angles, that have the greatest
%        range. 'lab' includes a string matrix with the angles for labeling.
%        'ang1' provides which angle has the greatest range.
%
%        meaning of the input arguments:
%          - direction_matrix: matrix including the elevation and azimuth
%                              angle of the different directions. The first
%                              row includes the elevation, the second row
%                              the azimuth. Each column represents a
%                              direction.
%          - n_directions: number of directions
%
%        For further information see: meshplot, meshing, contplot, cont_ing
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

% seperates angles
elev = direction_matrix(1,:);
azi = direction_matrix(2,:);

% evaluate minimum and maximum angles
minelev = min(elev);
maxelev = max(elev);
minazi = min(azi);
maxazi = max(azi);

% detect maximum range
elevdelta = maxelev - minelev;
azidelta = maxazi - minazi;

if (elevdelta<azidelta)
  ang1 = 2;
  ang2 = 1;
else 
  ang1 = 1;
  ang2 = 2;
end

% initialize vector for y-axis

ytic = direction_matrix(ang1,:);
pseudo_ytic = direction_matrix(ang2,:);

ytic = ytic';
pseudo_ytic = pseudo_ytic';

% sets the string matrix for labeling

label1 = str2mat(num2str(ytic(1)));
label2 = str2mat(' ');
label3 = str2mat(',');
label5 = str2mat(num2str(pseudo_ytic(1)));

for i=2:n_directions
  newlabel = str2mat(num2str(ytic(i)));
  label1 = str2mat(label1,newlabel);
  label2 = str2mat(label2,' ');
  label3 = str2mat(label3,',');
  newlabel = str2mat(num2str(pseudo_ytic(i)));
  label5 = str2mat(label5,newlabel);
end;

lab = [label1,label2,label3,label2,label5];

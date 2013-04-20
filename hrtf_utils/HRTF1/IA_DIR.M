function ia_dir(ia_level,ia_time,direction_matrix,n_directions)

% function ia_dir(ia_level,ia_time,direction_matrix,n_directions)
%
% IA_DIR labels the directions in a 2D plot of the interaural level difference
%        versus the interaural time difference.
%
%        meaning of the input arguments:
%        - ia_level: column vector including the interaural level difference
%                    for each direction
%        - ia_time: column vector including the interaural time difference
%                   for each direction
%        - direction_matrix: matrix including the elevation and azimuth
%                            angle of the different directions. The first
%                            row includes the elevation, the second row
%                            the azimuth. Each column represents a
%                            direction.
%        - n_directions: number of directions
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


[x,y]=ginput(1);
maxlevel = max(ia_level);
minlevel = min(ia_level);
yrange = maxlevel - minlevel;
ydelta  = yrange/20;

maxtime = max(ia_time);
mintime = min(ia_time);
xrange = maxtime - mintime;
xdelta = xrange/20;

j=0;
for i=1:n_directions
  xdiff=abs(ia_time(i)-x);
  ydiff=abs(ia_level(i)-y);
  if (xdiff<=xdelta)&(ydiff<=ydelta)
    j = j + 1;
    direct(j) = i;
    differences(j) = sqrt(xdiff^2+ydiff^2); 
  end;
end;

[which_diff,ind] = min(differences);
index = direct(ind);
text(ia_time(index),ia_level(index)+ydelta*1.5,['elev: ',num2str(direction_matrix(1,index))]);
text(ia_time(index),ia_level(index)+ydelta/2,['azi: ',num2str(direction_matrix(2,index))]);
hold on;
plot(ia_time(index),ia_level(index),'r+');
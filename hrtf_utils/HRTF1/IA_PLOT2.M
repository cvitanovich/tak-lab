function ia_plot2(level_difference,time_difference,direction_matrix,...
                  n_directions,frq_axis)

% function ia_plot2(level_difference,time_difference,direction_matrix,...
%                   n_directions,frq_axis)
%
% IA_PLOT2 plots a 3D plot of the interaural level difference versus the
%          interaural time difference for all evaluated frequencies.
%          The different directions are plotted in different colors. For
%          each frequency, the different directions are connected with lines,
%          too (like IA_PLOT1). Sliders are provided for changing the view
%          position of the 3D-plot.
%
%          meaning of the input arguments:
%          - level_difference: matrix including the interaural level
%                              differences. Each row represents one direction,
%                              each column represents the level differences
%                              for a certain frequency.
%          - time_difference: matrix including the interaural time differences.
%                             Each row represents one direction, each column
%                             represents the level differences for a certain
%                             frequency.
%          - direction_matrix: matrix including the elevation and azimuth
%                              angle of the different directions. The first
%                              row includes the elevation, the second row
%                              the azimuth. Each column represents a
%                              direction.
%          - n_directions: number of directions
%          - frq_axis: vector including the frequencies in [Hz], for which
%                      the interaural differences are evaluated.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter

colorarray = [n_directions];
for i=1:n_directions
  colorindex = rem(i,7);    % seven different colors
  if colorindex==1
     colorarray(i) = 'y';
  elseif colorindex==2
     colorarray(i) = 'm';
  elseif colorindex==3
     colorarray(i) = 'c';
  elseif colorindex==4
     colorarray(i) = 'r';
  elseif colorindex==5
     colorarray(i) = 'g';
  elseif colorindex==6
     colorarray(i) = 'b';
  elseif colorindex==0
     colorarray(i) = 'w';
  end;
end;

% plotting
n_iaural = length(frq_axis);

% draw line between different frequencies of one direction

h = figure;
for i=1:n_directions
  ia_level = level_difference(i,:);
  ia_time = time_difference(i,:);
  plot3(frq_axis,ia_time,ia_level,colorarray(i));
  hold on;
end;

% draw line between different directions and a single frequency

for i=1:n_iaural
  ia_level = level_difference(:,i);
  ia_time = time_difference(:,i);
  frq_point = frq_axis(i) * ones(1,n_directions);
  plot3(frq_point,ia_time,ia_level,'w');
  hold on;
end;

% labels lines with directions
% for first frequency

[min_time,min_pos] = min(time_difference(:,1));
[max_time,max_pos] = max(time_difference(:,1));

elev = direction_matrix(1,min_pos);
azi = direction_matrix(2,min_pos);

frq_point = frq_axis(1);
plot3(frq_point,min_time,level_difference(min_pos,1),'y+');
text('pos',[frq_point,min_time,level_difference(min_pos,1)],'color','y','string',['E: ',num2str(elev),' A: ',num2str(azi)]);


elev = direction_matrix(1,max_pos);
azi = direction_matrix(2,max_pos);

frq_point = frq_axis(1);
plot3(frq_point,max_time,level_difference(max_pos,1),'y+');
text('pos',[frq_point,max_time,level_difference(max_pos,1)],'color','y','string',['E: ',num2str(elev),' A: ',num2str(azi)]);

% labels for last frequency

[min_time,min_pos] = min(time_difference(:,n_iaural));
[max_time,max_pos] = max(time_difference(:,n_iaural));

elev = direction_matrix(1,min_pos);
azi = direction_matrix(2,min_pos);

frq_point = frq_axis(n_iaural);
plot3(frq_point,min_time,level_difference(min_pos,n_iaural),'y+');
text('pos',[frq_point,min_time,level_difference(min_pos,n_iaural)],'color','y','string',['E: ',num2str(elev),' A: ',num2str(azi)]);


elev = direction_matrix(1,max_pos);
azi = direction_matrix(2,max_pos);

frq_point = frq_axis(n_iaural);
plot3(frq_point,max_time,level_difference(max_pos,n_iaural),'y+');
text('pos',[frq_point,max_time,level_difference(max_pos,n_iaural)],'color','y','string',['E: ',num2str(elev),' A: ',num2str(azi)]);

xlabel('f / Hz');
ylabel('time difference / sec');
zlabel('level difference / dB');
title('interaural differences');

grid on;

axissetting = axis;
axissetting(1) = frq_axis(1);
axissetting(2) = frq_axis(n_iaural);
axis(axissetting);

v3dcntrl(h);
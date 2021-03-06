function sph_cont(level_contour, labels)

% function sph_cont(level_contour,labels)
%
% SPH_CONT plots contour lines on a sphere and labels them. 'level_contour'
%          is the matrix C returned by the function 'contourc'.
% labels: 0) none 1) draw, default
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum% (p) 09.09.1994 Joerg Schmatter

if (nargin <2) labels =1;
[rows,columns] = size(level_contour);
next_level = 1;
label_distance = 20;
while next_level<columns
  n_o_pairs = level_contour(2,next_level);
  first_pair = next_level + 1;
  last_pair = next_level + n_o_pairs;
  azi_line = level_contour(1,first_pair:last_pair);
  elev_line = level_contour(2,first_pair:last_pair);
  % change elevation and azimuth to rectangular coordinates
%  x_line = cos(elev_line*pi/180) .* sin(-azi_line*pi/180);
%  y_line = -cos(elev_line*pi/180) .* cos(azi_line*pi/180);
%  z_line = sin(elev_line*pi/180);
  [x_line,y_line,z_line]=sph2cart(azi_line/180*pi,elev_line/180*pi,ones(size(elev_line)));
 
  hold on;
  line('XData',x_line,'YData',y_line,'ZData',z_line,'Color','white');
  n_o_labels = fix(n_o_pairs/label_distance);
  text_delta = .05;
  if (labels ==1)
    for i=0:n_o_labels-1
      x_text = x_line(i*label_distance + 1);
      y_text = y_line(i*label_distance + 1);
      z_text = z_line(i*label_distance + 1);
      if x_text < 0
        x_text = x_text - text_delta;
      else
        x_text = x_text + text_delta;
      end;
      if y_text < 0
        y_text = y_text - text_delta;
      else
        y_text = y_text + text_delta;
      end;
      if z_text < 0
        z_text = z_text - text_delta;
      else
        z_text = z_text + text_delta;
      end;
      text('position',[x_text,y_text,z_text],'color','k','string',num2str(level_contour(1,next_level)));
    end;
    next_level = next_level + n_o_pairs + 1;

  end;
end;

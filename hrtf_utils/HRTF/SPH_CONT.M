function sph_cont(level_contour, labels)

% function sph_cont(level_contour)
%
% SPH_CONT plots contour lines on a sphere and labels them. 'level_contour'
%          is the matrix C returned by the function 'contourc'.
% labels: 0) none 1) draw (default)
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

if (nargin <2) labels =1; end;
[rows,columns] = size(level_contour);
next_level = 1;
label_distance = 40;			%was 20
while next_level<columns
  n_o_pairs = level_contour(2,next_level);
  first_pair = next_level + 1;
  last_pair = next_level + n_o_pairs;
  azi_line = level_contour(1,first_pair:last_pair);
  elev_line = level_contour(2,first_pair:last_pair);

  % change elevation and azimuth to rectangular coordinates

  x_line = cos(elev_line*pi/180) .* sin(-azi_line*pi/180 + pi/2);
  y_line = -cos(elev_line*pi/180) .* cos(azi_line*pi/180 + pi/2);
  z_line = sin(elev_line*pi/180 );
 
  hold on;
if (level_contour(1,next_level) == 0)
  line('XData',x_line,'YData',y_line,'ZData',z_line,'Color','white');
else
  line('XData',x_line,'YData',y_line,'ZData',z_line,'Color','cyan');
end
  if (labels ==1)
    n_o_labels = fix(n_o_pairs/label_distance);
    text_delta = .05;
    for i=0:2:n_o_labels-1
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
      text('position',[x_text,y_text,z_text],'color','y','string',num2str(level_contour(1,next_level)));
    end;
  end;
  next_level = next_level + n_o_pairs + 1;
end;

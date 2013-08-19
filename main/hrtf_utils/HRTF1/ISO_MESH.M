function iso_mesh(isofrequency,level_grid,level_color,level_contour,...
                  azi_interp,elev_interp)

% function iso_mesh(isofrequency,level_grid,level_color,level_contour,...
%                   azi_interp,elev_interp)
%
% ISO_MESH plots lines of the same level in a 3D-meshplot. The different
%          levels are displayed with a pseudocolor plot and additional
%          contour lines. The view position can be changed with sliders
%          controlling the azimuth and elevation angle.
%
%          Description of the input arguments:
%           - isofrequency: selected frequency that is plotted
%           - level_grid: matrix including the interpolated magnitudes in
%                         [dB]. Each column represents a certain azimuth angle,
%                         each row represents a certain elevation angle.
%           - level_color: matrix including the pseudocolors belonging to the
%                          interpolated magnitudes. Each column represents a
%                          certain azimuth angle, each row represents a 
%                          certain elevation angle.
%           - level_contour: matrix including the output arguments of the function
%                            'contourc'
%           - azi_interp: vector including the azimuth angles
%           - elev_interp: vector including the elevation angles
%
%           For further information see: iso_fit
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


mesh_h = figure;
mesh(azi_interp,elev_interp,level_grid,level_color);
grid on;

[rows,columns] = size(level_contour);
next_level = 1;
label_distance = 35;
while next_level<columns
  n_o_pairs = level_contour(2,next_level);
  first_pair = next_level + 1;
  last_pair = next_level + n_o_pairs;
  x_line = level_contour(1,first_pair:last_pair);
  y_line = level_contour(2,first_pair:last_pair);
  z_line = level_contour(1,next_level) * ones(1,n_o_pairs);
  hold on;
  line('XData',x_line,'YData',y_line,'ZData',z_line,'Color','white');
  n_o_labels = fix(n_o_pairs/label_distance);
  for i=0:n_o_labels-1
    x_text = level_contour(1,first_pair + i*label_distance);
    y_text = level_contour(2,first_pair + i*label_distance);
    z_text = level_contour(1,next_level) + 1;
    text(x_text,y_text,z_text,num2str(level_contour(1,next_level)));
  end;
  next_level = next_level + n_o_pairs + 1;
end;

xlabel('azimuth / deg');
ylabel('elevation / deg');
zlabel('level / dB');

axissetting = axis;
axissetting(1) = azi_interp(1);
axissetting(2) = azi_interp(length(azi_interp));
axissetting(3) = elev_interp(1);
axissetting(4) = elev_interp(length(elev_interp));
axis(axissetting);

%doesn't work 2/96 CHK:
%v3dcntrl(mesh_h);


function isoglobe(isofrequency,level_grid,level_color,level_contour,...
                  azi_interp,elev_interp)

% function isoglobe(isofrequency,level_grid,level_color,level_contour,...
%                   azi_interp,elev_interp)
%
% ISOGLOBE plots lines of the same level on a sphere. The different levels are
%          displayed with a pseudocolor plot and additional contour lines.
%          The sphere includes a net of the polar coordinates for a better
%          orientation. The lines of the same elevation and the lines of the
%          same azimuth have an angle resolution of 15°. The view position
%          can be changed with sliders controlling the azimuth and elevation
%          angle.
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
%           - elev_intep: vector including the elevation angles
%
%           For further information see: iso_fit
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


iso_radius = 0.98;

% changing directions to rectangular coordinates

disp('Transformation from sphere to rectangular coordinates');

iso_x = [];
iso_y = [];
iso_z = [];


for i=1:length(elev_interp)
  for j=1:length(azi_interp)
[iso_x(i,j),iso_y(i,j),iso_z(i,j)]=sph2cart(azi_interp(j)/180*pi,elev_interp(i)/180*pi,iso_radius);
  end;
end;


isoplot_h = figure;

% plotting of the approximated iso_points
Y_N = input('Plot surface (y/n) ','s');
if (Y_N == 'y')
  surf_h = surf(iso_x,iso_y,iso_z,level_color);
  set(surf_h,'edgecolor','none');
end;
% plotting of the sphere net
Y_N = input('Plot Long & Latt (y/n) ','s');
if (Y_N == 'y')
  polarnet;
end;

% plotting contour lines on the surface
% leave this out for now, CHK 2/96:
Y_N = input('Plot contours (y/n) ','s');
if (Y_N == 'y')
 sph_cont(level_contour);
end;

% sphere orientation

%text('pos',[1.5,0,0],'string','front','color','w');
text('pos',[0,-1.5,0],'string','right','color','w');
%text('pos',[-1.5,0,0],'string','rear','color','w');
text('pos',[0,1.5,0],'string','left','color','w');


% slides and view control
% this doesn't work CHK 2/96:
%v3dcntrl(isoplot_h);


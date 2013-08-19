function ia_globe(difference,direction_matrix,frq_axis,type)

% function ia_globe(difference,direction_matrix,frq_axis,type)
%
% IA_GLOBE plots a sphere with a combined plot of the iso - interaural level
%          or interaural time differences with additional contour lines and
%          pseudocolors for a single frequency. Sliders are provided to change
%          the view position of the sphere.
%
%          Description of the input arguments:
%          - difference: vector including the interaural time or level differences
%                        of the measured directions at a selected frequency
%          - direction_matrix: matrix including the elevation and azimuth
%                              angle of the different directions. The first
%                              row includes the elevation, the second row
%                              the azimuth. Each column represents a
%                              direction.
%          - frq_axis: vector including the frequencies in [Hz], for which
%                      the interaural differences are evaluated.
%          - type: separates between interaural 'level_difference' and
%                  'time_difference' for plotting.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


[n_directions,n_iaural] = size(difference);

if n_iaural~=1
  [selected_frq,selected_line] = ia_s_frq(frq_axis,n_iaural);
  plot_difference = difference(:,selected_line);
else
  selected_frq = frq_axis(1);
  plot_difference = difference;
end;

% generate higher resolution of the sampled data array at the selected frequency

[level_grid,level_color,level_contour,azi_interp,...
 elev_interp] = ia_fit(plot_difference,direction_matrix,n_directions,type);

isoglobe(selected_frq,level_grid,level_color,level_contour,...
         azi_interp,elev_interp);

if strcmp(type,'level_difference')
  title(['interaural level difference in [dB]        ',num2str(selected_frq),' Hz']);
elseif strcmp(type,'time_difference');
  title(['interaural time difference in [10^(-6) sec]        ',num2str(selected_frq),' Hz']);
end;
axis off
prcntrl
view(-90,90);
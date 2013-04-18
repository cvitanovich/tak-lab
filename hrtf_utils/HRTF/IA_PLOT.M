function ia_plot(level_difference,time_difference,direction_matrix,...
                 frq_axis,single_plane)

% function ia_plot(level_difference,time_difference,direction_matrix,...
%                  frq_axis,single_plane)
%
% IA_PLOT plots for a horizontal plane the interaural time difference and the
%         interaural level difference about the changing azimuth angle for a
%         single frequency.
%
%         Description of the input arguments:
%         - level_difference: column vector or matrix including the interaural
%                             level differences for the different directions
%         - time_difference: column vector or matrix including the interaural
%                            level differences for the different directions
%         - direction_matrix: matrix including the elevation and azimuth
%                             angle of the different directions. The first
%                             row includes the elevation, the second row
%                             the azimuth. Each column represents a
%                             direction.
%         - frq_axis: including the frequency in [Hz], for which the
%                     interaural differences are evaluated
%         - single_plane: including the information if the interaural differences
%                         are evaluated in a single horizontal plane or not
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter



[n_directions,n_iaural] = size(level_difference);

if n_iaural~=1
  [selected_frq,selected_line] = ia_s_frq(frq_axis,n_iaural);

  ia_level = level_difference(:,selected_line);
  ia_time = time_difference(:,selected_line);

else
  ia_time = time_difference;
  ia_level = level_difference;
  selected_frq = frq_axis(1);
end;

if strcmp(single_plane,'elev')
  figure;
  subplot(2,1,1);
  plot(direction_matrix(2,:),ia_time,'y-',direction_matrix(2,:),ia_time,'y+');
  axissetting = axis;
  axissetting(1) = direction_matrix(2,1);
  axissetting(2) = direction_matrix(2,n_directions);
  axis(axissetting);
  title(['elevation: ',num2str(direction_matrix(1,1)),'°     frequency: ',num2str(selected_frq),' Hz']);
  ylabel('interaural time difference in [sec]');
  xlabel('azimuth in deg');

  subplot(2,1,2);
  plot(direction_matrix(2,:),ia_level,'y-',direction_matrix(2,:),ia_level,'y+');
  axissetting = axis;
  axissetting(1) = direction_matrix(2,1);
  axissetting(2) = direction_matrix(2,n_directions);
  axis(axissetting);
  title(['elevation: ',num2str(direction_matrix(1,1)),'°     frequency: ',num2str(selected_frq),' Hz']);
  ylabel('interaural level difference in [dB]');
  xlabel('azimuth in deg');
else
  disp('function only available for a single horizontal plane');
end;


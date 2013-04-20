function ia_plot1(level_difference,time_difference,direction_matrix,...
                  n_directions,frq_axis,single_plane)

% function ia_plot1(level_difference,time_difference,direction_matrix,...
%                   n_directions,frq_axis,single_plane)
%
% IA_PLOT1 plots the interaural level difference versus the interaural time
%          difference. If the interaural differences are evaluated in a
%          single horizontal plane the evaluated points are connected with
%          a line, else they are plotted without any connection. A menu offers
%          the possibility to label the plotted points with the elevation and
%          azimuth belonging to them.
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
%         - n_directions: number of different directions
%         - frq_axis: including the frequency in [Hz], for which the
%                     interaural differences are evaluated
%         - single_plane: including the information if the interaural differences
%                         are evaluated in a single horizontal plane or not
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter




n_iaural = length(frq_axis);

if n_iaural~=1
  [selected_frq,selected_line] = ia_s_frq(frq_axis,n_iaural);

  ia_level = level_difference(:,selected_line);
  ia_time = time_difference(:,selected_line);

else             % only single frequency available
  ia_level = level_difference;
  ia_time = time_difference;
  selected_frq = frq_axis(1);
end;

% plotting

figure;

if strcmp(single_plane,'elev')
  plot(ia_time,ia_level,'y+',ia_time,ia_level,'y-');
else
  plot(ia_time,ia_level,'y+');
end;

xlabel('interaural time difference / sec');
ylabel('interaural level difference / dB');
title(['frequency: ',num2str(selected_frq),' Hz']);

% direction labeling

label_choice = 0;
while label_choice~=2
  label_choice=menu('Direction labeling',...
                    'Select',...
                    'Close');
  if label_choice==1
    ia_dir(ia_level,ia_time,direction_matrix,n_directions);
  end;
end;
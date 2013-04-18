% IA_CNTRL provides a menu for different kind of plots of interaural level
%          differences versus interaural time differences.
%
%          - Single frequency: 2D plot of level difference versus time
%                              difference, if all directions are situated in
%                              a single horizontal plane.
%          - All frequencies:  3D plot of level difference versus time
%                              difference for all evaluated frequencies.
%
%          For further information see: ia_plot1, ia_plot2
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


plot_choice = 0;

while (plot_choice~=1)&(plot_choice~=2)
  plot_choice = menu(' ',...
                     'Single frequency',...
                     'All frequencies');
  if plot_choice==1
    if strcmp(single_plane,'elev')
      ia_plot1(level_difference,time_difference,direction_matrix,...
               n_directions,frq_axis,single_plane);
    else
      disp('Function only available for datas of a single horizontal plane');
    end;
    return;
  elseif plot_choice==2
    ia_plot2(level_difference,time_difference,direction_matrix,...
             n_directions,frq_axis);
    return;
  end;
end;
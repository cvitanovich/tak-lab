% ISOLYSE guides through the visualisation methods for isolines of
%         transfer functions of the external ears.
%         The following methods are available:
%             - Load data:  loading a catalog of transfer functions
%             - Sphere plot: Combined 3D plot of the isolevels in a 3D sphere
%                            with contour lines and pseudocolors
%             - 3D-Plane plot: Combined 3D plot of the isolevels above a plane
%                              in three dimensions with contour lines and
%                              pseudocolors
%             - Exit: exits the menu
%
%         The loaded data must have the filetype 3 (real- and imaginary
%         part) and stored in the MTLSIG - Format.
%         The measured directions are linear interpolated with the function
%         'iso_fit'. If the file includes two channels the channel has to
%         be selected after choosing the plotting method. Each plot displays
%         several isolines for one selected frequency.
%
%         For further information see: MTLSIG-manuals, iso_fit, isoglobe,
%                                      iso_mesh, isoplot
%           
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter



clc;
choice = 0;
%dataloaded = 0; 

while choice~=5
  choice=menu('Visualizing isolines',...
              'Load data',...
              'Sphere plot Level',...
              'Sphere plot Phase',...
              '3D-Plane plot',...
              'Exit');
  if (choice==1)
    [filename,pathname] = uigetfile('*.mtl','Reading MTLSIG-Format');
    if (filename~=0)
      clc;  				% clear workspace
      filename = [pathname filename];
      [file_format,filetype,first_channel,second_channel,first_phase,...
       second_phase,direction_matrix,time_delay,n_directions,n_channels,...
       n_lines,sampling_rate,begin_of_signal,end_of_signal,...
       xaxis,comment] = read(filename);
      if (file_format=='MTLSIG')&(filetype==3)
        dataloaded = 1;
      else
        disp('File hasn''t MTLSIG-Format or incorrect filetype!');
        dataloaded = 0;
      end;
    end;
  elseif (choice==2)
    if dataloaded==1
      if n_channels==2
         selected_channel = 0;
         while selected_channel==0
           selected_channel = menu('Select channel',...
                                   'First channel',...
                                   'Second channel',...
				   'Difference');
         end;
      else
        selected_channel = 1;
      end;
      plot_kind = 'sp_level';
      plot_kind
      isoplot(filename,abs2db(first_channel),abs2db(second_channel),selected_channel,...
              direction_matrix,n_directions,n_lines,sampling_rate,...
              begin_of_signal,end_of_signal,plot_kind);
    else
      disp('no data loaded...'); 
    end;

 elseif (choice==3)
    if dataloaded==1
      if n_channels==2
         selected_channel = 0;
         while selected_channel==0
           selected_channel = menu('Select channel',...
                                   'First channel',...
                                   'Second channel',...
				   'Difference');
         end;
      else
        selected_channel = 1;
      end;
      plot_kind = 'sp_phase';
      
      isoplot(filename,angle(first_channel)/(pi),angle(second_channel)/(pi),...
              selected_channel,...
              direction_matrix,n_directions,n_lines,sampling_rate,...
              begin_of_signal,end_of_signal,plot_kind);
    else
      disp('no data loaded...'); 
    end;

  elseif (choice==4)
    if dataloaded==1
      if n_channels==2
         selected_channel = 0;
         while selected_channel==0
           selected_channel = menu('Select channel',...
                                   'First channel',...
                                   'Second channel');
         end;
      else
        selected_channel = 1;
      end;
      plot_kind = 'plane';
      isoplot(filename,abs2db(first_channel),abs2db(second_channel),selected_channel,...
              direction_matrix,n_directions,n_lines,sampling_rate,...
              begin_of_signal,end_of_signal,plot_kind);
    else
      disp('no data loaded...'); 
    end;
  end;
end;


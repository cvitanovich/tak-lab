% IA_LYSE guides through the different visualisation methods of interaural
%         level differences and interaural time differences.
%         The following methods are available:
%         - Load data: loading a catalog of transfer functions or interaural
%                      differences stored with the MTLSIG - format.
%         - Evaluation: evaluate interaural level difference and interaural 
%                       time difference for a single frequency or all loaded
%                       frequencies in selected steps.
%         - Save data: save evaluated interaural differences in the MTLSIG -
%                      format
%         - Sphere-Plot of level differences: plot interaural level difference
%                                             on a sphere
%         - Sphere-Plot of time differences: plot interaural level difference
%                                            on a sphere
%         - Plot level and time differences: plot for a horizontal plane the
%                                            interaural level and time difference
%                                            for a single frequency
%         - Level versus Time differences: plot the interaural level difference
%                                          versus the interaural time difference
%                                          for a single frequency in a 2D plot
%                                          or for all frequencies in a 3D plot.
%         - Exit: exits the menu
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter

clc;
choice = 0;
data_loaded = 0;
evaluated = 0;

while choice~=8
  choice=menu('Interaural analyse','Load data',...
              'Evaluation',...
              'Save data',...
              'Sphere-Plot of level difference',...
              'Sphere-Plot of time difference',...
              'Plot level and time difference',...
              'Level versus Time difference',...
              'Exit');

 if (choice==1)
    [filename,pathname] = uigetfile('*.mtl','Reading MTLSIG-Format');
    if (filename~=0)
      clc;  % clear workspace
      filename = [pathname filename];
     [file_format,filetype,first_channel,second_channel,first_phase,...
       second_phase,direction_matrix,time_delay,n_directions,n_channels,...
       n_lines,sampling_rate,begin_of_signal,end_of_signal,...
       xaxis,comment] = read(filename);
      if (file_format~='MTLSIG')
        disp('File has not MTLSIG-Format !');
        dataloaded = 0;
      elseif (filetype~=3)&(filetype~=11)
        disp('File has incorrect filetype !');
        dataloaded = 0;
      elseif (filetype==11)
        level_difference = first_channel;
        time_difference = first_phase;
        frq_axis = xaxis;
        clear first_channel;
        clear second_channel;
        clear xaxis;
        evaluated = 1;
        dataloaded = 1;
        single_plane = planedct(direction_matrix);
      elseif (filetype==3)
        if n_channels==2
          evaluated = 0;
          dataloaded = 1;
          single_plane = planedct(direction_matrix);
        else
          disp(' ');
          disp('For interaural analysis the file must include two channels !');
          disp('Load new file !');
          dataloaded = 0;
        end;
      end;
    end;
  elseif (choice==2)          % evaluation
    if (filetype==3)&(dataloaded==1)
       frq_choice = 0;

       while (frq_choice~=1)&(frq_choice~=2)
         frq_choice = menu(' ',...
                           'Single frequency',...
                           'All frequencies');
         if frq_choice==1
           [level_difference,time_difference,frq_axis] = ia_eval1(first_channel,...
                                                                  second_channel,n_directions,...
                                                                  n_lines,sampling_rate,begin_of_signal,...
                                                                   end_of_signal);

         elseif frq_choice==2
         [level_difference,time_difference,frq_axis] = ia_evalx(first_channel,...
                                                                  second_channel,n_directions,...
                                                                  n_lines,sampling_rate,begin_of_signal,...
                                                                  end_of_signal,xaxis);

         end;
       end;

      evaluated = 1;
    end;
  elseif (choice==3)          % save data
    if evaluated==1
      save_name = [];
      while isempty(save_name)
        save_name = input('Filename for saving: ','s');
      end;
      save_name = lower(save_name);
      ia_save(save_name,level_difference,time_difference,direction_matrix,...
              n_directions,sampling_rate,frq_axis,comment);
    end;
  elseif (choice==4)            % plot interaural level difference
    if evaluated==1
      if strcmp(single_plane,'elev')
        disp('Only single horizontal plane available; choose another file');
      else
        ia_globe(level_difference,direction_matrix,frq_axis,'level_difference');
      end;
    else
      disp('Choose ''Evaluation'' at first !');
    end;
  elseif (choice==5)            % plot interaural time difference
    if evaluated==1
      if strcmp(single_plane,'elev')
        disp('Only single horizontal plane available; choose another file');
      else
        ia_globe(time_difference,direction_matrix,frq_axis,'time_difference');
      end;
    else
      disp('Choose ''Evaluation'' at first !');
    end;
  elseif (choice==6)            % plot interaural level and time difference for a single plane
    if evaluated==1
      ia_plot(level_difference,time_difference,direction_matrix,frq_axis,single_plane);
    end;
  elseif (choice==7)
    if evaluated==1
        if length(frq_axis)==1
          ia_plot1(level_difference,time_difference,direction_matrix,...
                   n_directions,frq_axis,single_plane);
        else
          ia_cntrl;
        end;
    end;
  end;
end;
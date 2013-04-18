% ANALYSE guides through the different visualisation methods of impulse
%         responses and transfer functions of the external ears.
%         The following methods are available:
%           - Load data: loading a catalog of impulse responses or
%                        transfer functions with the MTLSIG or SIGMAP -format
%           - Magnitude: displaying the magnitude over the time, the
%                        frequency or a psycho-acoustic scale
%           - Phase: displaying the phase of the transfer functions over
%                    the frequency or a psycho-acoustic scale
%           - Phase delay / Group delay: displaying the the phase or group
%                                        delay of signals with the filetype
%                                        5 or 6.
%           - Mesh: 3D mesh displaying of the impulse responses or transfer
%                   functions of all directions
%           - Pseudocolor / Contour: pseudocolor and contour plot of the
%                                    transfer functions of all directions
%           - Contour: contour plot of the transfer functions of all
%                      directions
%           - Movie: 2D or 3D movie of the impulse responses or transfer
%                    functions
%           - Exit: exits the menu
%
%         After loading a catalog of signals, the duration or bandwidth can
%         be reduced and for frequency signals the kind of scaling (linear
%         or logarithm) of the frequency axis must be choosen.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


choice = 0;
dataloaded = 0;

while choice~=9
  choice=menu('Visualizing Tools',...
              'Load data',...
              'Magnitude',...
              'Phase',...
              'Phase delay / Group delay',...
              'Mesh',...
              'PseudoColor - Contour',...
              'Contour',...
              'Movie',...
              'Exit');
  if (choice==1)			% loading file
     [filename,pathname] = uigetfile('*.mtl;*.sc*;*.sr*','READING data');
      if (filename~=0)
        clc;  % clear workspace
        filename = [pathname filename];
        [file_format,filetype,first_channel,second_channel,first_phase,...
         second_phase,direction_matrix,time_delay,n_directions,n_channels,...
         n_lines,sampling_rate,begin_of_signal,end_of_signal,...
         xaxis,comment] = read(filename);
        dataloaded = 1;
        [first_channel,second_channel,first_phase,second_phase,n_lines,...
         begin_of_signal,end_of_signal,xaxis] = cut_sig(filetype,...
         first_channel,second_channel,first_phase,second_phase,n_directions,...
         n_channels,n_lines,sampling_rate,begin_of_signal,end_of_signal,...
         xaxis);
        if ( (filetype>=3)&(filetype<=6) )    % linear frequency signals
          scaling = frqscale;
        else
          scaling = 'lin';
        end;
      end;
  elseif (choice==2)			% magnitude plots
    if dataloaded==1
      if scaling=='lin'
        magnplot(file_format,filetype,filename,first_channel,second_channel,...
                 direction_matrix,n_directions,n_channels,n_lines,xaxis);
      elseif scaling=='log'
        magn_log(file_format,filetype,filename,first_channel,second_channel,...
                 direction_matrix,n_directions,n_channels,n_lines,xaxis);
      end;
    else disp('no data loaded...'); 
    end;
  elseif (choice==3)			% phase plots
    if dataloaded==1
      if scaling=='lin'
        ph_plot(file_format,filetype,filename,first_channel,second_channel,...
                first_phase,second_phase,direction_matrix,n_directions,...
                n_channels,n_lines,xaxis);
      elseif scaling=='log'
        ph_log(file_format,filetype,filename,first_channel,second_channel,...
               first_phase,second_phase,direction_matrix,n_directions,...
               n_channels,n_lines,xaxis);
      end;
    else disp('no data loaded...'); 
    end;   
  elseif (choice==4)			% phase delay or group delay plots
            % only available for filetype 5 (phase delay) and 6 (group delay)
            % in MTLSIG format
    if dataloaded==1
      if scaling=='lin'
        delay(filetype,filename,first_phase,second_phase,direction_matrix,...
              n_directions,n_channels,n_lines,xaxis);
      else
        delaylog(filetype,filename,first_phase,second_phase,direction_matrix,...
                 n_directions,n_channels,n_lines,xaxis);
      end;
    else disp('no data loaded...');
    end;
  elseif (choice==5)			% meshplots
    if dataloaded==1
      meshplot(file_format,filetype,filename,first_channel,second_channel,...
               first_phase,second_phase,direction_matrix,n_directions,...
               n_channels,n_lines,xaxis,scaling);
    else disp('no data loaded...'); 
    end;
  elseif (choice==6)			% pseudocolor contours
    if dataloaded==1
     psc_cont(file_format,filetype,filename,first_channel,second_channel,first_phase,...
              second_phase,direction_matrix,n_directions,n_channels,...
              n_lines,xaxis,scaling);
    else disp('no data loaded...'); 
    end;
  elseif (choice==7)			% contour plots
    if dataloaded==1
     contplot(file_format,filetype,filename,first_channel,second_channel,first_phase,...
              second_phase,direction_matrix,n_directions,n_channels,...
              n_lines,xaxis,scaling);
    else disp('no data loaded...');
  end;
  elseif (choice==8)			% movie
    if dataloaded==1
      if file_format=='MTLSIG'
        mo_cntrl;
      else
        disp('file must include MTLSIG format');
      end;
    else disp('no data loaded...'); 
    end;
  end;
end;


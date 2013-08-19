function ph_log(file_format,filetype,filename,first_channel,second_channel,...
                first_phase,second_phase,direction_matrix,n_directions,...
                n_channels,n_lines,xaxis)

% function ph_log(file_format,filetype,filename,first_channel,second_channel,...
%                 first_phase,second_phase,direction_matrix,n_directions,...
%                 n_channels,n_lines,xaxis)
%
% PH_LOG Plots the phase, phase delay or group delay of the signals
%        over 'xaxis' with a logarithm x-axis scaling. Signals with the
%        filetype 4,5 and 6 plot 'f_phase' and 's_phase'. Signals with the
%        filetype 3,7,8,9 and 10 plot the evaluated phase of 'first_channel'
%        and 'second_channel'.
%
%        meaning of the input arguments:
%         - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                        labeling
%         - filetype: number of identification for the signal type
%         - filename: name of the file
%         - first_channel: matrix with the absolute magnitude of the
%                          signals first channel
%         - second_channel: matrix with the absolute magnitude of the signals
%                           second channel. If the file includes only one
%                           channel, this matrix is empty.
%         - first_phase: phase, phase delay or group delay of the first
%                        channel
%         - second_phase: phase, phase delay or group delay of the second
%                         channel, if the file includes two channels. Else
%                         the variable is empty.
%         - direction_matrix: matrix including the elevation and azimuth
%                             angle of the different directions. The first
%                             row includes the elevation, the second row
%                             the azimuth. Each column represents a
%                             direction.
%         - n_directions: number of directions
%         - n_channels: number of channels
%         - n_lines: number of lines
%         - xaxis: vector including the points of sampling
%
%        For further informations see: MTLSIG manuals
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

if n_directions>=4
  max_plot_in_figure = 4;         % maximum number of subplots in each figure
else
  max_plot_in_figure = n_directions;
end;

if filetype==2    % time signal
  disp('time signals haven''t any phase'); 
elseif ( (filetype==3)|(filetype==7)|(filetype==8)|(filetype==9)|(filetype==10) )
  if n_channels==1          % windows for the phase
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = 180/pi*unwrap(angle(first_channel(signal_to_plot,:)));
        semilogx(xaxis,signal);
        p_layout(file_format,filetype,filename,k,'transfer function: ',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  else
    % plot first channel;   windows for the phase
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = 180/pi*unwrap(angle(first_channel(signal_to_plot,:)));
        semilogx(xaxis,signal);
        p_layout(file_format,filetype,filename,k,'transfer function first channel: ',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end

   % plot second channel;         windows for the phase
   % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = 180/pi*unwrap(angle(second_channel(signal_to_plot,:)));
        semilogx(xaxis,signal);
        p_layout(file_format,filetype,filename,k,'transfer function second channel: ',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
   end
  end % endif number_of_channels
elseif filetype==4        % signals with magnitude and phase description
  if n_channels==1          % windows for the phase
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = first_phase(signal_to_plot,:);
        semilogx(xaxis,signal);
        p_layout(file_format,filetype,filename,k,'transfer function: ',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  else
    % plot first channel;   windows for the phase
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = first_phase(signal_to_plot,:);
        semilogx(xaxis,signal);
        p_layout(file_format,filetype,filename,k,'transfer function first channel: ',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end

   % plot second channel;         windows for the phase
   % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = second_phase(signal_to_plot,:);
        semilogx(xaxis,signal);
        p_layout(file_format,filetype,filename,k,'transfer function second channel: ',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  end % endif number_of_channels
else
  disp('signals include phase delay or group delay, use other function'); 
end % endif signaltype seperation
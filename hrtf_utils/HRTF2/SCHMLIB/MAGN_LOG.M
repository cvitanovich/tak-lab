function magn_log(file_format,filetype,filename,first_channel,second_channel,...
                  direction_matrix,n_directions,n_channels,n_lines,xaxis)

% function magn_log(file_format,filetype,filename,first_channel,second_channel,...
%                   direction_matrix,n_directions,n_channels,n_lines,xaxis)
%
% MAGN_LOG Plots the signals magnitude (first_channel,second_channel)
%          over 'xaxis' with a logarithm x-axis.
%
%          meaning of the input arguments:
%          - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                         labeling
%          - filetype: number of identification for the signal type
%          - filename: name of the file
%          - first_channel: matrix with the absolute magnitude of the signals
%                           first channel
%          - second_channel: matrix with the absolute magnitude of the signals
%                            second channel. If the file includes only one
%                            channel, this matrix is empty.
%          - direction_matrix: matrix including the elevation and azimuth angle
%                              of the different directions. The first row
%                              includes the elevation, the second row the
%                              azimuth. Each column represents a direction.
%          - n_directions: number of directions
%          - n_channels: number of channels
%          - n_lines: number of lines
%          - xaxis: vector including the points of sampling
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
  if n_channels==1
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = first_channel(signal_to_plot,:);
        semilogx(xaxis,signal);
        t_layout(file_format,filetype,filename,k,'impulse response',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  else                                               % windows for the first channel
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = first_channel(signal_to_plot,:);
        semilogx(xaxis,signal);
        t_layout(file_format,filetype,filename,k,'impulse response: first channel',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
                                                     % windows for the second channels
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = second_channel(signal_to_plot,:);
        semilogx(xaxis,signal);
        t_layout(file_format,filetype,filename,k,'impulse response: second channel',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  end
else                                            % all frequency signals
  if n_channels==1       % windows for the magnitude
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    all_signals = abs2db(first_channel);
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = all_signals(signal_to_plot,:);
        semilogx(xaxis,signal);
        m_layout(file_format,filetype,filename,k,'transfer function: magnitude',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  else
    % plot first channel
    % controling the number of figures and the number of subplots for one figure
    [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    all_signals = abs2db(first_channel);
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = all_signals(signal_to_plot,:);
        semilogx(xaxis,signal);
        m_layout(file_format,filetype,filename,k,'transfer function first channel: magnitude',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
   
   % plot second channel
   % controling the number of figures and the number of subplots for one figure
   [number_of_figures,figure_list,number_of_subplots] = fg_cntrl(max_plot_in_figure,n_directions,xaxis);
    signal_to_plot = 0;
    all_signals = abs2db(second_channel);
    for j=1:number_of_figures
      figure(figure_list(j));
      n_o_s = number_of_subplots(j);
      for k=1:n_o_s
        subplot(max_plot_in_figure,1,k);
        signal_to_plot = signal_to_plot + 1;
        signal = all_signals(signal_to_plot,:);
        semilogx(xaxis,signal);
        m_layout(file_format,filetype,filename,k,'transfer function second channel: magnitude',signal_to_plot,direction_matrix,xaxis,n_lines);
      end
      zoom_sub(xaxis);
    end
  end % end number_of_channels
end % end time signal / frequency signal
function delay(filetype,filename,first_phase,second_phase,...
               direction_matrix,n_directions,n_channels,n_lines,xaxis)

% function delay(filetype,filename,first_phase,second_phase,...
%                direction_matrix,n_directions,n_channels,n_lines,xaxis)
%
% DELAY plots the phase delay or group delay of signals with the filetype 5 and 6.
%       The x-axis has a linear scaling.
%
%       - filetype: number of identification for the signal type
%       - filename: name of the file
%       - first_phase: phase or group delay of the first channel
%       - second_phase: phase or group delay of the second channel, if the
%                       file includes two channels. Else the variable is
%                       empty.
%       - direction_matrix: matrix including the elevation and azimuth angle
%                           of the different directions. The first row
%                           includes the elevation, the second row the azimuth.
%                           Each column represents a direction.
%       - n_directions: number of directions 
%       - n_channels: number of channels
%       - n_lines: number of lines
%       - xaxis: vector including the points of sampling
%
%       For further information see MTLSIG - manuals.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


if n_directions>=4
  max_plot_in_figure = 4;         % maximum number of subplots in each figure
else
  max_plot_in_figure = n_directions;
end;

if (filetype~=5)&(filetype~=6)
  disp('loaded filetype doesn''t include phase delay or group delay');
else
 if n_channels==1       % windows for the magnitude
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
        plot(xaxis,signal);
        if (filetype==5)
          delayout(filetype,filename,k,'phase delay  ',signal_to_plot,direction_matrix,xaxis,n_lines);
        else
          delayout(filetype,filename,k,'group delay  ',signal_to_plot,direction_matrix,xaxis,n_lines);
        end;
      end
      zoom_sub(xaxis);
    end
 else
    % plot first channel
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
        plot(xaxis,signal);
        if (filetype==5)
          delayout(filetype,filename,k,'phase delay first_channel  ',signal_to_plot,direction_matrix,xaxis,n_lines);
        else
          delayout(filetype,filename,k,'group delay first channel  ',signal_to_plot,direction_matrix,xaxis,n_lines);
        end;
      end
      zoom_sub(xaxis);
    end
   
   % plot second channel
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
        plot(xaxis,signal);
        if (filetype==5)
          delayout(filetype,filename,k,'phase delay second channel  ',signal_to_plot,direction_matrix,xaxis,n_lines);
        else
          delayout(filetype,filename,k,'group delay second channel  ',signal_to_plot,direction_matrix,xaxis,n_lines);
        end; 
      end;
      zoom_sub(xaxis);
    end;
 end; % end number_of_channels
end;

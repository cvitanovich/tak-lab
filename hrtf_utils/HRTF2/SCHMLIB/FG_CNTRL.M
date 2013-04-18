function [number_of_figures,figure_list,number_of_subplots]...
         = fg_cntrl(max_plot_in_figure,n_directions,xaxis)

% function [number_of_figures,figure_list,number_of_subplots]...
%         = fg_cntrl(max_plot_in_figure,n_directions,xaxis)
%
% FG_CNTRL evaluates the necessary number of figures for plotting and opens
%          them. The function returns the number of figures, a list of
%          handels to each figure and the number of subplots in each figure.
%          The size of the subplots depends only on 'max_plot_in_figure'.
%          So it is possible that figures are not used totaly.
%          'n_directions' includes the number of directions and 'xaxis'
%          is a vector including the points of sampling.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


number_of_figures = ceil(n_directions / max_plot_in_figure);
n_o_s = n_directions;
for i=1:number_of_figures
  figure_list(i) = figure;
  n_o_s = n_o_s - max_plot_in_figure;
  if n_o_s<0
    number_of_subplots(i) = n_o_s + max_plot_in_figure;
  else
    number_of_subplots(i) = max_plot_in_figure;
  end
end

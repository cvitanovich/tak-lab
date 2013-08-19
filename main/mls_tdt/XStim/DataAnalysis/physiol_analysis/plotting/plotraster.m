function [] = plotraster(hrasterfig,...
   spikedata,...
   xdimval,...
   ydimval,...
   repnum,...
	tot_dur,...
   stim_dur,...
   numreps,...
   xoffset,...
   yoffset,...
   row_space,...
   col_space,...
   marker_color);

if nargin < 13 | isempty(marker_color)
    marker_color = 'c';
end

figure(hrasterfig)
ind_under = find(spikedata <= stim_dur);
ind_over = find(spikedata>stim_dur & spikedata<tot_dur);
spikes_under = spikedata(ind_under);
spikes_over = spikedata(ind_over);
 
 [x y] = time2graph(spikes_under,...
    xdimval(ind_under),...
    ydimval(ind_under),...
    repnum, tot_dur, numreps, xoffset, yoffset, row_space, col_space);
ph = plot(x,y, '.','color',marker_color);
set(ph, 'MarkerSize', 3);
if(~isempty(spikes_over));
   [x y] = time2graph(spikes_over,...
      xdimval(ind_over),...
      ydimval(ind_over),...
      repnum,...
      tot_dur,...
      numreps,...
      xoffset,...
      yoffset,...
      row_space,...
      col_space);
   ph = plot(x,y, 'r.');
   set(ph, 'MarkerSize', 3);
end

return;

function [x, y] = time2graph(spikedata, xdimval, ydimval, repnum, tot_dur, numreps, ...
                            xoffset, yoffset, row_space, col_space);
% TIME2GRAPH  converts spike times into raster grid plot values
%
%  time2graph(data_array, dur, nreps, ...
%                          xoffset, yoffset, row_space, col_space);
%  data_array has the format [spike_time rep xval_num, yval_num]
%
% converts from spike time, rep number, x-dim value number, y-dim value number into
% graph coordinates.  time, rep, xval_num, and yval_num can be col vectors
% This function is written for the case where all spikes in a raster grid are plotted
% on a single graph.
% see also draw_grid()

originx = (xdimval-1)*(tot_dur+2*xoffset+col_space) + xoffset;
originy = (ydimval-1)*(numreps+2*yoffset+row_space) + yoffset;
x = spikedata+originx;
y = (repnum-1)+originy;

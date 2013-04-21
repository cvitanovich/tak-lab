function [hrastfig,hrastaxes] = makerasterfig(position,...
   xvals,...
   yvals,...
   xtext,...
   ytext,...
   xoffset,...
   yoffset,...
   numreps,...
   tot_dur,...
   stim_dur,...
   startcount,...
   endcount,...
   row_space,...
   col_space)

background = 'black';
foreground = 'blue';
textcolor = 'yellow';

hrastfig = figure('Position',position,...
   'Color',background,...
   'Name','Raster Plot');
hrastaxes = axes;
set(hrastaxes,'Position',[0.10 0.05 .95 0.95]);
hold on;
axis off;

max_ygraphval = (length(yvals))*(numreps+2*yoffset+row_space); 
max_xgraphval = (length(xvals)-1)*(tot_dur+2*xoffset+col_space);

for y = 1:length(yvals)
  for x = 1:length(xvals)
     xpos = (x-1)*(tot_dur+2*xoffset+col_space);
     ypos = (y-1)*(numreps+2*yoffset+row_space); 
     draw_rast(xpos, ypos, xoffset, yoffset, numreps, tot_dur, stim_dur, startcount, endcount, foreground);
     if (x==1) 
        th = text(xpos-5, ypos+(numreps+yoffset)/2, num2str(yvals(y)));
        set(th, 'HorizontalAlignment', 'right');
        set(th, 'FontSize', 8, 'Color', textcolor);
     end;
     if (y == 1 & length(xvals)>1)
        th = text(xpos+(tot_dur+xoffset)/2, y-5, num2str(xvals(x))); 
        set(th, 'HorizontalAlignment', 'center');
        set(th, 'VerticalAlignment', 'top');
        set(th, 'FontSize', 8, 'Color', textcolor);
     end;
  end;
end;


return;

function x = draw_rast(xpos, ypos, xoffset, yoffset, numreps, dur, stim_dur, startcount, endcount, foreground)
% draws the borders of a raster at the specified location


bracket_length = 50;

line1x = [xpos xpos];
line1y = [ypos ypos+(numreps+2*yoffset-1)];

line2x = [xpos+dur+2*xoffset-1 xpos+dur+2*xoffset-1];
line2y = line1y;

bracket1x = [xpos xpos+bracket_length];
bracket1y = [ypos ypos];

bracket2x = bracket1x;
bracket2y = (ypos+numreps+2*yoffset-1)*[1 1];

bracket3x = [xpos+dur+2*xoffset-1-bracket_length xpos+dur+2*xoffset-1];
bracket3y = bracket1y;

bracket4x = bracket3x;
bracket4y = bracket2y;

stim_offx = (xpos+stim_dur+xoffset)*[1 1];
stim_offy = line1y; 

if (~isempty(startcount))
   stcntx = (xpos+startcount+xoffset)*[1 1];  % start count x coordinate
   stcnty = line1y;
end;

if (~isempty(endcount))
   endcntx = (xpos+endcount+xoffset)*[1 1];
   endcnty = line1y;
end;

new_width = .1;
l = line(line1x, line1y);
set(l, 'linewidth', new_width, 'Color', foreground);
l = line(line2x, line2y);
set(l, 'linewidth', new_width, 'Color', foreground);
l = line([line1x(1) line2x(1)], [line1y(1) line1y(1)]);
set(l, 'linewidth', new_width, 'Color', foreground);
%l = line([line1x(1) line2x(1)], [line2y(1) line2y(1)]);
%set(l, 'linewidth', new_width, 'Color', foreground);
l = line(bracket1x, bracket1y);
set(l, 'linewidth', new_width, 'Color', foreground);
l = line(bracket2x, bracket2y);
set(l, 'linewidth', new_width, 'Color', foreground);
l = line(bracket3x, bracket3y);
set(l, 'linewidth', new_width, 'Color', foreground);
l = line(bracket4x, bracket4y);
set(l, 'linewidth', new_width, 'Color', foreground);
if ~(dur==stim_dur)   % rev cor tests have no raster overrun
   l = line(stim_offx, stim_offy);
   set(l, 'linewidth', new_width, 'Color', foreground);
end;

if (~isempty(startcount))
   l = line(stcntx, stcnty);
   set(l, 'linewidth', new_width);
   set(l, 'color', 'green');
end;

if (~isempty(endcount))
   l = line(endcntx, endcnty);
   set(l, 'linewidth', new_width);
   set(l, 'color', foreground);
end;

return;

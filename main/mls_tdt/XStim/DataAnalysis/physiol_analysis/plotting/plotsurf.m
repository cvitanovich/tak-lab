function [disp_az, disp_el, disp_data] = plotsurf(xaxis, yaxis, data);
% PLOTSURF plots 2-d data using pcolor and corrects axes
%   [disp_az, disp_el, disp_data] = plotsurf(xaxis, yaxis, spacedata);
% NOTE: THIS VERSION HAS BEEN UPDATED TO SHIFT NUMERICAL VALUES OF AXES
%       RATHER THAN JUST THE LOCATION OF THE TICKS

% display parts not shown using pcolor and adjust axes
y_step = min(diff(yaxis));
x_step = min(diff(xaxis));
n_xaxis = size(data,2);
n_yaxis = size(data,1);
if (size(xaxis,1) > 1) xaxis = xaxis'; end;
if (size(yaxis,1) > 1) yaxis = yaxis'; end;

if (y_step>0)
   disp_y = [yaxis yaxis(n_yaxis) + y_step];
else
	disp_y = [yaxis(1) - y_step yaxis];
end;

if (x_step>0)
	disp_x =  [xaxis  xaxis(n_xaxis)  + x_step];
else
	disp_x = [xaxis(1) - x_step xaxis];
end;

disp_data = [data data(:,n_xaxis)];
disp_data = [disp_data; disp_data(n_yaxis,:)];
	
% realign axes
yaxis_shift = abs(y_step/2);
xaxis_shift = abs(x_step/2);
	
pcolor(disp_x-xaxis_shift,disp_y-yaxis_shift,disp_data);

shading flat;

return;

% realign axes  (cut out because its better to just shift the data as above)
yaxis_shift = abs(y_step/2);
xaxis_shift = abs(x_step/2);
ca = gca;
xtl = get(ca,'xticklabel');
ytl = get(ca,'yticklabel');
set(ca,'xtick',get(ca,'xtick')+xaxis_shift);
set(ca,'ytick',get(ca,'ytick')+yaxis_shift);
set(ca,'xticklabel',xtl);
set(ca,'yticklabel',ytl);	

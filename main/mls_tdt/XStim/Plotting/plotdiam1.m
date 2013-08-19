function [disp_az, disp_el, disp_data] = plotdiam1(dir,data);
% PLOTDIAM1 plots space data using pcolor and corrects axes
%   [] = plotdiam1(dir, spacedata);

contournum = 0; 

% display parts not shown in using pcolor and adjust axes
elstep = 5;
azstep = 5;
disp_el = [-90:5:90];
disp_az = [-90:5:90];
n_azi = size(disp_el,2);
n_ele = size(disp_az,2);


warning('off')
disp_data = ones(n_azi, n_ele) ./ nan;
warning('on')

for i = 1:size(data,2)
   row = find(disp_el == dir(1,i));
   col = find(disp_az == dir(2,i));
   disp_data(row,col) = data(i);
end

% shift to realign axes
disp_el = disp_el - elstep/2;
disp_az = disp_az - azstep/2;
% add row & column for pcolor
disp_data = [disp_data disp_data(:,n_azi)];
disp_data = [disp_data; disp_data(n_ele,:)];
disp_az = [disp_az disp_az(n_azi)];
disp_el = [disp_el disp_el(n_ele)];

if (contournum)
   contour(disp_az, disp_el, disp_data, contournum, contourcolor);
else
   pcolor(disp_az,disp_el,disp_data);
end;

shading flat;
ca = gca;
set(ca,'dataaspectratio', [1 1 1])
set(ca,'Color','black');

% set axis limits and shade unavailable HRTFs
axis([-90 90 -90 90])	
patch([-90 -90 0], 	[0 90 90], 	[0.9 0.9 0.9]);
patch([-90 -90 0], 	[0 -90 -90], [0.9 0.9 0.9]);
patch([90 90 0], [0 90 90], 	[0.9 0.9 0.9]);
patch([90 90 0], [0 -90 -90], 	[0.9 0.9 0.9]);

return;
	
% realign axes
yaxis_shift = elstep/2;
xaxis_shift = azstep/2;
ca = gca;
xtl = get(ca,'xticklabel');
ytl = get(ca,'yticklabel');
set(ca,'xtick',get(ca,'xtick')+xaxis_shift);
set(ca,'ytick',get(ca,'ytick')+yaxis_shift);
set(ca,'xticklabel',xtl);
set(ca,'yticklabel',ytl);
set(ca,'dataaspectratio', [1 1 1])
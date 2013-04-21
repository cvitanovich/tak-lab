function [disp_az, disp_el, disp_data] = plotdiam(azi, ele, spacedata, contournum, contourcolor);
% PLOTDIAM plots space data using pcolor and corrects axes
%   [] = plotdiam(azi, ele, spacedata);

if (nargin<4) 
  contournum = 0; 
  contourcolor = 'w';
end;

% display parts not shown in using pcolor and adjust axes
elstep = min(diff(ele));
if isempty(elstep) 	elstep = 5;	end
azstep = min(diff(azi));
if isempty(azstep)	azstep = 5;	end
n_azi = size(spacedata,2);
n_ele = size(spacedata,1);
disp_el = [ele ele(n_ele) + elstep];
disp_az = [azi  azi(n_azi)  + azstep];

% shift to realign axes
disp_el = disp_el - elstep/2;
disp_az = disp_az - azstep/2;

disp_data = [spacedata spacedata(:,n_azi)];
disp_data = [disp_data; disp_data(n_ele,:)];

figure
if (contournum)
   contour(disp_az, disp_el, disp_data, contournum, contourcolor);
else
   pcolor(disp_az,disp_el,disp_data);
end;

shading flat;
ca = gca;
set(ca,'dataaspectratio', [1 1 1])

% set axis limits and shade unavailable HRTFs
axis([-90 90 -90 90])	
patch([-90 -90 0], 	[0 90 90], 		[0.9 0.9 0.9]);
patch([-90 -90 0], 	[0 -90 -90], 	[0.9 0.9 0.9]);
patch([90 90 0], 		[0 90 90], 		[0.9 0.9 0.9]);
patch([90 90 0], 		[0 -90 -90], 	[0.9 0.9 0.9]);

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
	

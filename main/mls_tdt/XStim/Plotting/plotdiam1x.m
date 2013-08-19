function [disp_az, disp_el, disp_data] = plotdiam1(azi, ele, spacedata, interpflag, contournum, contourcolor);
% PLOTDIAM1 plots space data using pcolor and corrects axes
%   [] = plotdiam1(azi, ele, spacedata, interpflag);

if nargin < 4
    interpflag = 1;
end
if (nargin<5) 
  contournum = 0; 
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

%*************************
if interpflag
for i = 1:length(ele)-1
   for j = 1:length(azi)-1
      if isnan(disp_data(i,j))
         if i ==1
            if j == 1
               temp = [disp_data(i+1,j),disp_data(i,j+1)];
            else
               temp = [disp_data(i+1,j),disp_data(i,j-1),disp_data(i,j+1)];
            end
         else
            if j == 1
               temp = [disp_data(i-1,j),disp_data(i+1,j),disp_data(i,j+1)];
				else
               temp = [disp_data(i-1,j),disp_data(i+1,j),disp_data(i,j-1),disp_data(i,j+1)];
            end
         end
         
         ind = find(isnan(temp) ==0);
         if length(ind) >=3
            disp_data(i,j) = mean(temp(ind));
         end
      end
   end
end
end


%**********************************
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
	

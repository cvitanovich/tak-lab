function flatplot(v, MAP, contourF, colour, res, step, az, el, slider, titlestr, nofigure,minv,maxv)

% FLATPLOT interpolated three dimensional colour, mesh or contour plot  
%
% FLATPLOT(v, MAP, contourF, colour, res, step, az, el, slider, titlestr, nofigure,minv,maxv)
% plot transfer function on sphere
% v:
% MAP: 
% line : number of line to plot
% contourF: 0 no contour lines,1 with contour lines
% colour: 0 no colour 1 with colour
% res: spatial resolution (5 deg)
% step: stepwidth for contourlines
% az: azimuth of viewing angle
% el: elevation of viewing angle
% slider: with(1) or without(0) slider default (1)
% titlestring: string for toplin of figure
% nofigure: if 1 flatplot does not create a new figure, slider are disabled
% minv: minimum value for colour table
% maxv: maximum value for colour table

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 26.10.1995

% Calls         Toolbox
%
%v3dcntrl
%prcntrl
%(sphereinterp)


if (nargin < 3) contourF=1; end;
if (nargin < 4) colour=1; end;
if (nargin < 5) res=5; end;
if (nargin < 6) step=0; end;
if (nargin < 7) az=-90; el=90; end;
if (nargin < 8) el=90; end;
if (nargin < 9) slider=1; end;
if (nargin < 10) titlestr=''; end;
if (nargin < 11) nofigure=0; end;
if (nargin < 12) minv=min(min(v));end;
if (nargin < 13) maxv=max(max(v));end;

if (nofigure==0) plot_handle=figure; end;

fig_title=titlestr;
azi_step=res;
ele_step=res;
dir_matrix= MAP;
v_min=min(min(v));
v_max=max(max(v));
if (step ~=0) c_step=step;    
else  c_step=1; end;      
min_c=c_step*fix(v_min/c_step);
max_c=c_step*fix(v_max/c_step);
n_c=(max_c-min_c)/c_step+1;
cc=linspace(min_c',max_c,n_c);

% interpolation

%disp('Interpolation');

min_azi = min(dir_matrix(2,:));
max_azi = max(dir_matrix(2,:));
delta=abs(max_azi-min_azi);
if (delta>180) delta=abs(delta-360); end;

% added by K. Hartung 19.2.95
%if (delta<=45) min_azi=-max_azi; end; 
% added k. Hartung 19.2.95if (delta<=45) max_azi=min_azi+360; end;
% added k. Hartung 19.2.95

% find number of different elevations
ele=dir_matrix(1,:);
ele=sort(ele);
d=[1 diff(ele)];
ele=ele(find(d));
min_ele = min(ele);
max_ele = max(ele);
n_ele=length(ele);   % find number of different azimuths

azi_hires=min_azi:azi_step:max_azi;
ele_hires=min_ele:ele_step:max_ele;

for i=1:n_ele  
 azi_index=find(ele(i)==dir_matrix(1,:)); 
 n_azi=length(azi_index); 
 azi=dir_matrix(2,azi_index);  
 [azi,index]=sort(azi);  
 azi_index=azi_index(index);   
 v_lores=v(azi_index);   
 v_c=zeros(1,n_azi+2);  
 azi_c=zeros(1,n_azi+2); 
 azi_c(1)=azi(n_azi)-360; 
 azi_c(2:n_azi+1)=azi; 
 azi_c(n_azi+2)=azi(1)+360;  
 azi=azi_c;  
 v_c(1)=v_lores(n_azi); 
 v_c(2:n_azi+1)=v_lores;
 v_c(n_azi+2)=v_lores(1);
 v_hires(i,:)= (interp1(azi,v_c,azi_hires))';
end; 

v_hires=interp2(azi_hires,ele,v_hires,azi_hires,ele_hires); 

% surface colors
if (colour)
 %disp('Surface colors')
 surf_h=surface(azi_hires,ele_hires,v_hires);
 level_color = get(surf_h,'Cdata');
end

% contour lines

if (contourF)
 %disp('Evaluate contour lines');
 v_contour = contourc(azi_hires,ele_hires,v_hires,cc);
end;

% plot
%disp ('Plotting');

colormap('jet');

if (colour)   
  surface_handle=surf(azi_hires,ele_hires,v_hires,level_color);
  set(surface_handle,'edgecolor','none') 
else
%  surface_handle=surf(azi_hires,ele_hires,v_hires,(v_hires.*0));  
  set(surface_handle,'edgecolor','none');  
  set(surface_handle,'facecolor','k');
end;

caxis([minv,maxv]);
axis([min(azi_hires),max(azi_hires),min(ele_hires),max(ele_hires),....
min(min(v_hires)),max(max(v_hires))]);
hold on
%caxis([minv,maxv]);
shading interp

% plotting of the sphere net grid on the sphere

%polarnet2([min(azi),max(azi)],[min(ele),max(ele)],2,1);
% plotting contour lines on the surface
%if (contour) sph_cont(v_contour); end;

if (contourF) 
 cs=contour3(azi_hires,ele_hires,v_hires); 
 clabel(cs);
end;

hold off
ylabel('Elevation [0]');xlabel('Azimuth [0]');
zlabel('Level [dB]');title(titlestr);
% sphere orientation
%text('pos',[1.5,0,0],'string','front','color','w');
%text('pos',[0,-1.5,0],'string','right','color','w');
%text('pos',[-1.5,0,0],'string','rear','color','w');
%text('pos',[0,1.5,0],'string','left','color','w');
% slides and view control
if ((slider==1) & (nofigure==0)) 
     v3dcntrl(plot_handle);
  else    
  colorbar('vert'); 
end;

az=0;el=90;
view(az,el);
title(fig_title);

%switch on printer control
if (nofigure~=1)
  prcntrl; 
end




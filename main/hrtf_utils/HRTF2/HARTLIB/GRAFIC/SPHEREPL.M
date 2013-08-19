function SPHEREPL(v, MAP, contour, colour, res, step, az, el, slider, nose, minv, maxv, fig_title, bar_title)

% SPHEREPL interpolated three dimensional colour and/or contour plot on sphere 
% SPHEREPL(v, MAP, contour, colour, res, step, az, el, slider, nose, minv, maxv,fig_title,bar_title)
%
% plot interaural transfer function on sphere
% filename
% line : number of line to plot
% contour: 0 no contour lines,1 with contour lines
% colour: 0 no colour 1 with colour
% res: spatial resolution (5 deg)
% step: stepwidth for contourlines
% az: azimuth of viewing angle
% el: elevation of viewing angle
% slider: with(1) or without(0) slider default (1)
% nose: (1) noseplot on 
% minv: minimum value for colour table
% maxv: maximum value for colour table

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 26.10.1995


if (nargin < 3) contour=1; end;
if (nargin < 4) colour=1; end;
if (nargin < 5) res=5; end;
if (nargin < 6) step=0; end;
if (nargin < 7) az=+90; el=0; end;
if (nargin < 8) el=0; end;
if (nargin < 9) slider=1; end;
if (nargin <10) nose=0; end;
if (nargin <11) maxv=max(max(v)); end;
if (nargin <12) minv=min(min(v)); end;
if (nargin <13) fig_title=''; end;
if (nargin <14) bar_title=''; end;
plot_handle=figure;

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

disp('Interpolation');

min_azi = min(dir_matrix(2,:));
max_azi = max(dir_matrix(2,:));

delta=abs(max_azi-min_azi);
if (delta>180) delta=abs(delta-360); end;% added by K. Hartung 19.2.95
%if (delta<=45) min_azi=-max_azi; end; % added k. Hartung 19.2.95
if (delta<=45) max_azi=min_azi+360; end;% added k. Hartung 19.2.95

% find number of different elevations
ele=dir_matrix(1,:);
ele=sort(ele);
d=[1 diff(ele)];
ele=ele(find(d));

min_ele = min(ele);
max_ele = max(ele);

n_ele=length(ele);

% find number of different azimuths

azi_hires=min_azi:azi_step:max_azi;
ele_hires=min_ele:ele_step:max_ele;

for k=1:n_ele
  azi_index=find(ele(k)==dir_matrix(1,:));
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

  v_hires(k,:)= (interp1(azi,v_c,azi_hires))';

end;

v_hires2=interp2(azi_hires,ele',v_hires,azi_hires,ele_hires');

% surface colors
if (colour)
disp('Surface colors')

%size(azi_hires)
%azi_hires;
%size(ele_hires)
%ele_hires;
%size(v_hires2)
%v_hires;
%surf_h=surface(azi_hires,ele_hires,(v_hires2-minv)/(maxv-minv));
surf_h=surface(azi_hires,ele_hires,v_hires2);
caxis([minv,maxv]);

level_color = get(surf_h,'Cdata');

%max(max(level_color))
%min(min(level_color))
%pause

end;

% contour lines
if (contour)
  disp('Evaluate contour lines');
  v_contour = contourc(azi_hires,ele_hires,v_hires2,cc);
end;

% transformation to cartesian coordinates

iso_radius=0.98;

n_azi_hires=length(azi_hires);
n_ele_hires=length(ele_hires);

x=zeros(n_ele_hires, n_azi_hires);
y=zeros(n_ele_hires, n_azi_hires);
z=zeros(n_ele_hires, n_azi_hires);

r=zeros(1,n_azi_hires)+iso_radius;
clc

disp('Transformation from sphere to cartesian')

if (nose==1)
  m_v_hires=mean(mean(v_hires));
  min_v_hires=min(min(v_hires));
  v_hires=v_hires-min_v_hires*zeros(size(v_hires));
  max_v_hires=max(max(v_hires));
  v_hires=v_hires/max_v_hires;
end;

for i=1:n_ele_hires
  ele_i=zeros(1,n_azi_hires)+ele_hires(i);
  if (nose==1)
    [x(i,:),y(i,:),z(i,:)]=sph2cart(azi_hires/180*pi,ele_i/180*pi,v_hires2(i,:)*0.1+1.0);
  else 
    [x(i,:),y(i,:),z(i,:)]=sph2cart(azi_hires/180*pi,ele_i/180*pi,r);
  end;
end;

% plot

disp ('Plotting');

if (colour) 
  surface_handle=surf(x,y,z,level_color);
  set(surface_handle,'edgecolor','none')
  axis off;
 else
 % surface_handle=surf(x,y,z,(z.*0));
 % set(surface_handle,'edgecolor','none')
 % set(surface_handle,'facecolor','k')
end;
caxis([minv,maxv]);

shading interp

% plotting of the sphere net grid on the sphere
% used to be polarnet2 - whatever that is CHK 10/17/96
%polarnet([min(azi), max(azi)], [min(ele), max(ele)], 2, 1);
polarnet;

if (slider)
     view(30,30); 
     v3dcntrl(plot_handle);
  else 
     h=colorbar('vert');
     t=get(h,'title');
     set(t,'string',bar_title);
     set(t,'FontSize',12);
     set(gca,'Xtick',[]);
     set(gca,'Ytick',[]);
     axis off;
end;

% plotting contour lines on the surface
% 'zero' indicates no labels on contours

if (contour) sph_cont(v_contour,0); end;

% sphere orientation

%get(gca,'FontSize')
set(gca,'FontSize',16);
%get(gca,'FontSize')

set(gca,'Xtick',[-1 0 1]);
set(gca,'Ytick',[-1 0 1]);
set(gca,'Ztick',[-1 0 1]);

set(gca,'Xlim',[-1 1]);
set(gca,'Ylim',[-1 1]);
set(gca,'Zlim',[-1 1]);

set(gca,'TickDir','in');

set(gca,'XtickLabels',[]);
set(gca,'YtickLabels',[]);
set(gca,'ZtickLabels',[]);

if (az>135) p=1; else p=-1; end; 	
%text('pos',[1.0,0,p*1.2],'string','front','color','w','FontSize',16,'HorizontalAlignment','center');
%text('pos',[0,-1.0,-p*1.2],'string','right','color','w','FontSize',16,'HorizontalAlignment','center');
%text('pos',[-1.0,0,-p*1.2],'string','rear','color','w','FontSize',16,'HorizontalAlignment','center');
%text('pos',[0,1.0,p*1.2],'string','left','color','w','FontSize',16,'HorizontalAlignment','center');

colormap('hot')

%set(gcf,'PaperOrientation','landscape');
%set(gcf,'PaperUnits','centimeters');
%set(gcf,'PaperType','a4letter');
%p = [2.5    2.5   26   16];
%set(gcf,'PaperPosition',p);

%hidden on
% slides and view control
 
set(gca,'AspectRatio',[1,1])

view(az,el)

%drawnow('discard');

%axis on
%grid on

%switch on printer control
%prcntrl

end;

%title(fig_title);
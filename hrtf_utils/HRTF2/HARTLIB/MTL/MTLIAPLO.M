function mtliaplot(filename, line, level, contour, colour, res, step, az, el)

%function mtliaplot(filename, line, level, contour, colour, res, step, az, el)
% plot interaural transfer function on sphere
% filename
% line : number of line to plot
% level: 0 (time difference) 1 level difference
% contour: 0 no contour lines,1 with contour lines
% colour: 0 no colour 1 with colour
% res: spatial resolution (5 deg)
% step: stepwidth for contourlines
% az: azimuth of viewing angle
% el: elevation of viewing angle

if (nargin < 3) level=0; end;
if (nargin < 4) contour=1; end;
if (nargin < 4) contour=1; end;
if (nargin < 5) colour=1; end;
if (nargin < 6) res=5; end;
if (nargin < 7) step=0; end;

if (nargin < 8) az=-90; el=90; end;
if (nargin < 9) el=90; end;

plot_handle=figure;

azi_step=res;
ele_step=res;

disp ('Load data');
% read header

[filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
           comment1, comment2] = mtlrh(filename);

% read frequency information

frq_axis = mtlrfrq(filename);
frq=frq_axis(line);

% read dir_matrix

dir_matrix= mtlrdir(filename,n_directions);

% read one line for all frequencies

ild=zeros(1,n_directions);
itd=zeros(1,n_directions);


for i=1:n_directions

 [ild(i),itd(i)]=mtlrch(filename, i, line);
end;
if (level==0) v=itd*10^(6);
 else v=ild; end;

v_min=min(v);
v_max=max(v);


if(level==0 ) 

	if (step ~=0) c_step=step; 
   		else  
         c_step=100; end;  %in microseconds
        
elseif (step ~=0) c_step=step;
       else  
         c_step=1;  %in dB
     end;
end;


min_c=c_step*fix(v_min/c_step);
max_c=c_step*fix(v_max/c_step);
n_c=(max_c-min_c)/c_step+1;
cc=linspace(min_c,max_c,n_c);

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
 

end;



% surface colors
if (colour)
disp('Surface colors')

surf_h=surface(azi_hires,ele_hires,v_hires);
level_color = get(surf_h,'Cdata');
end;

% contour lines
if (contour)
disp('Evaluate contour lines');

v_contour = contourc(azi_hires,ele_hires,v_hires,cc);

end;

% transformation to cartesian coordinates

iso_radius=0.98;

n_azi_hires=length(azi_hires);
n_ele_hires=length(ele_hires);

x=zeros(n_ele_hires, n_azi_hires);
y=zeros(n_ele_hires, n_azi_hires);
z=zeros(n_ele_hires, n_azi_hires);

r=zeros(1,n_azi_hires)+iso_radius;

disp('Transformation from sphere to cartesian')
for i=1:n_ele_hires
ele_i=zeros(1,n_azi_hires)+ele_hires(i);
  [x(i,:),y(i,:),z(i,:)]=sph2cart(azi_hires/180*pi,ele_i/180*pi,r);
end;

% plot

disp ('Plotting');

if (colour) 
  surface_handle=surf(x,y,z,level_color);
  set(surface_handle,'edgecolor','none')
 else
  surface_handle=surf(x,y,z,(z.*0));
  set(surface_handle,'edgecolor','none')
  set(surface_handle,'facecolor','k')
end;

% plotting of the sphere net grid on the sphere
polarnet2([min(azi),max(azi)],[min(ele),max(ele)],2,1);

% plotting contour lines on the surface

if (contour) sph_cont(v_contour); end;

% sphere orientation

text('pos',[1.5,0,0],'string','front','color','w');
text('pos',[0,-1.5,0],'string','right','color','w');
text('pos',[-1.5,0,0],'string','rear','color','w');
text('pos',[0,1.5,0],'string','left','color','w');


% slides and view control
v3dcntrl(plot_handle);
shading interp
view(az,el)

if (level==1) title(['interaural level difference in [dB]        ',num2str(frq),' Hz']);
else title(['interaural time difference in [10^(-6) sec]        ',num2str(frq),' Hz']);
end;
axis off

%switch on printer control
prcntrl


end;

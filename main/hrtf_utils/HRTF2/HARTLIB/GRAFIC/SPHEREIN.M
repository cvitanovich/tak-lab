function [v_hires,MAP_hires]=spherein(v,MAP,res)

% SPHEREIN interpolation on a sphere
%
% [v_hires,MAP_hires]=SPHEREINTERP(v,MAP,res)
%
% v: values to be interpolated
% MAP: coordinates of values
% res: resolution
% v_hires: interpolated values
% MAP_hires: coordinates of interplated values

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 16.09.1996 derived from flatplot and a routine of J. Schmatter
 
azi_step=res;
ele_step=res;

dir_matrix= MAP;
v_min=min(min(v));
v_max=max(max(v));




% interpolation%disp('Interpolation');


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
n_ele=length(ele);
% find number of different azimuths
azi_hires=min_azi:azi_step:max_azi;
ele_hires=min_ele:ele_step:max_ele;
for i=1:n_ele  azi_index=find(ele(i)==dir_matrix(1,:)); 
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
MAP_hires=meshgrid(azi_hires,ele_hires);

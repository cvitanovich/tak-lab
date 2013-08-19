function [azi, ele, data] = array2diamond(v, MAP, forplot);

% ARRAY2DIAMOND stuffs spatial information from 1D array to diamond matrix
%
%               [azi, ele, data] = array2diamond(v, MAP, forplot)
%
%               Reads the coordinates matrix and uses this to place
%               the array values in their correct spatial locations.
%               This works for double-polar coordinates that
%               are *evenly spaced* but not necessarily in any order.
%               Output locations not available under double polar
%               coordinates are set to NaN.
%               azi and ele are axes for the data grid 
%               the coordinate values.
%               MAP = a 2 by n array of coordinates [el; az]
%               v   = is a n by 1 array of values at each spatial
%                     locations (e.g. ILD values at a given freq)
%               forplot = 1 if output data should be suitable for
%                         plotting with pcolor (adds extra row and column and
%                         recalibrates axes accordingly) else 0

%swrap =0;
%colour =1;

if (nargin<3) forplot = 0; end;

dir_matrix= MAP;

% stuff data in rectangular matrix
%disp('Stuffing data');

min_azi = min(dir_matrix(2,:));
max_azi = max(dir_matrix(2,:));
min_ele = min(dir_matrix(1,:));
max_ele = max(dir_matrix(1,:));


% find number of different elevations and find increment
ele=dir_matrix(1,:);
ele=sort(ele);
d=[1 diff(ele)];
ele=ele(find(d));
n_ele=length(ele);

% set el increment
ele_step = max(d);

% find number of different azimuths and find az increment
% ASSUMES THIS IS DOUBLE POLAR DATA SO ALL AZIMUTHS ARE EQUALLY SPACED!!!!!!
azi   = dir_matrix(2,:);
azi   = sort(azi);
d     = [1 diff(azi)];
azi   = azi(find(d));
n_azi = length(azi);

azi_step = max(d);

% form data matrix
data = NaN*zeros(n_ele,n_azi);

% fill data matrix
m = length(v);
for i = 1:m
  row = (dir_matrix(1,i) - min_ele + azi_step)/azi_step;
  col = (dir_matrix(2,i) - min_azi + ele_step)/ele_step;
  data(row,col) = v(i);
end


% enlarge everything by one row & one column so the entire area is plotted
% this is needed when using pcolor for plotting results
% not needed for image
if (forplot)
	ele = [ele ele(n_ele) + azi_step];
	azi = [azi azi(n_azi) + ele_step];
	data = [data data(:,n_azi)];
	data = [data; data(n_ele,:)];
end;

return; 



























































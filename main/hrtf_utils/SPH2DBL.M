function [dir1, dir2] = sph2dbl(p1, p2)

%SPH2DBL.m	converts spherical coordinates into double polar coordinates
%		[dir1, dir2] = sph2dbl(p1, p2):  (in degrees)
%		input arguments: saz, sel 
%		output arguments: daz, del
%	or:	[dir1] = sph2dbl(p1):  (in degrees)
%		input argument: direction matrix (2 x n_dirs; el then az)
%		output argument: direction matrix (2 x n_dirs; del then daz)

if (nargin == 1)
  el = p1(1,:).*pi/180;
  az = p1(2,:).*pi/180;
%  disp('SPH2DBL: calculating double polar direction matrix');
elseif (nargin ==2)
  el = p1*pi/180;
  az = p2*pi/180;
%  disp('SPHDBL: calculating double polar azimuth and elevation');
else
  error('SPH2DBL: incorrect number of input arguments');
  return;
end;

% convert to cartesian coordinates:
x = cos(el) .* sin(az);
y = sin(el);
z = cos(el) .* cos(az);

%convert to double polar coordinates:
del = (el);
daz = (asin(x));

[value,index] = find(daz > pi/2);
daz(index) = pi - daz(index);

[value,index] = find(daz < -pi/2);
daz(index) = -pi - daz(index);

if (nargin == 1)
  dir1 = zeros(size(p1));
  dir1(1,:) = round(del*180/pi);
  dir1(2,:) = round(daz*180/pi);
else
  dir1 = round(del*180/pi);
  dir2 = round(daz*180/pi);
end;

%disp (['el : ' num2str(p1) ' del : ' num2str(dir1)]);
%disp (['az : ' num2str(p2) ' daz : ' num2str(dir2)]);

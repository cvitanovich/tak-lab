function [dir1, dir2] = dbl2sph(d1, d2)

%DBL2SPH.m	converts double polar into spherical coordinates
%		[dir1, dir2] = dbl2sph(d1, d2):  (in degrees)
%		input arguments: daz, del 
%		output arguments: saz, sel
%	or:	[dir1] = dbl2sph(p1):  (in degrees)
%		input argument: direction matrix (2 x n_dirs; del then daz)
%		output argument: direction matrix (2 x n_dirs; el then az)


% convert to radians
DegtoRad = pi/180;
RadtoDeg = 180/pi;

if (nargin == 1)
  del = d1(1,:) * DegtoRad;
  daz = d1(2,:) * DegtoRad;
%  disp('DBL2SPH: calculating double polar direction matrix');
elseif (nargin ==2)
  del = d1 * DegtoRad;
  daz = d2 * DegtoRad;
%  disp('DBL2SPH: calculating spherical azimuth and elevation');
else
  error('DBL2SPH: incorrect number of input arguments');
  return;
end;

% convert to cartesian
sel = del;
saz = asin(sin(daz)./cos(del));

if (saz == -1) saz = -1.57079; end;
if (saz == 1)  saz =  1.57079; end;

% back to degrees:
sa = round(saz * RadtoDeg);
se = round(sel * RadtoDeg);

if (nargin == 1)
  dir1 = zeros(size(d1));
  dir1(1,:) = se;
  dir1(2,:) = sa;
else
  dir1 = se;
  dir2 = sa;
end;

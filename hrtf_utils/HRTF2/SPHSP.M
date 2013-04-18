
% SPHSP Spherical Sample Points
% maps a set of point evenly space in elev and angle to another
% set of point which is evenly space in distance along the sphere.

el_start = -pi/2; el_end = pi/2;
el_incr = pi/30
az_start = 0; az_end = pi;
az_incr = pi/30;

% generate a set of points evenly spaced in radians elev and azim
el = el_start:el_incr:el_end;
az = az_start:az_incr:az_end;
[EL, AZ] = meshgrid(el, az);

% transform each point to another point in polar coordinates that will
% be evenly spaced.

new_az = asin(cos(AZ)./cos(EL));

real_az = new_az(isreal(new_az));
real_EL = EL(isreal(new_az));

hemisphere(real_az, real_EL); colormap(jet);
hidden off;
graph_rot

[X, Y, Z] = sph2cart(real(new_az), EL, ones(size(EL)));
figure;
surf(X,Y,Z); colormap(jet);

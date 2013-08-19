% POLARNET plots a net of lines, like a globe on an already existing sphere,
%          for better spatial orientation. POLARNET is based on the function
%          'sphere', but uses a higher resolution, for better drawing approximation.
%          The lines have an angle distance of 15° in azimuth and elevation
%          direction, and are plotted in black.
%
%          For further information see: isoplot
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter
% plot iso-elevation-planen = 180;              
%sphere resolution
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.

theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi);
cosphi(1) = 0; 
cosphi(n+1) = 0;
sintheta = sin(theta); 
sintheta(1) = 0; 
sintheta(n+1) = 0;
all_x = cosphi*cos(theta);
all_y = cosphi*sintheta;
all_z = sin(phi)*ones(1,n+1);
x = [];
y = [];
z = [];
for i=1:5:n+1   
  x = [x;all_x(i,:)];
  y = [y;all_y(i,:)];
  z = [z;all_z(i,:)];
end;
hold on;
plot3(x',y',z','w.');
% plot iso-azimuth-planen = 96;
              %sphere resolution
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi);
cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta);
sintheta(1) = 0;
sintheta(n+1) = 0;
all_x = cosphi*cos(theta);
all_y = cosphi*sintheta;
all_z = sin(phi)*ones(1,n+1);
u = [];
v = [];
w = [];
for i=1:4:n+1 
  u = [u;all_x(:,i)];
  v = [v;all_y(:,i)];
  w = [w;all_z(:,i)];
end;
hold on;
plot3(u',v',w','k.');
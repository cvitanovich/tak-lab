function stretch(f)
% STRETCH blows up the current axis by factor f.
%         STRETCH(f)

a=axis;
dx=f*abs(a(2)-a(1))/2;
dy=f*abs(a(4)-a(3))/2;
a(1)=a(1)-dx;
a(2)=a(2)+dx;
a(3)=a(3)-dy;
a(4)=a(4)+dy;
axis(a);

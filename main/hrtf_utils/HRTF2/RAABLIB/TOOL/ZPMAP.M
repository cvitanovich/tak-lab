function zpmap(b,a,gtitle,mode);
%ZPMAP  ZPMAP(b,a,gtitle,mode) plots the poles and
%       zeros of the transfer function H(z) specified
%       by b and a.
%
%               b(1)+ b(2)*z^(-1) + ...
%       H(z) = ---------------------------
%               a(1)+ a(2)*a^(-1) + ...
%
% Note: if any of the values in a or b are complex,
% then they are plotted direct!

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 15.08.1994 A. Raab

if nargin<4 mode='c'; end;
if nargin<3 gtitle=''; end;
if mode(1)=='f' pc='.r'; zc='.g'; else pc='xr'; zc='og'; end;

if ~(any(imag(a)~=0) | any(imag(b)~=0))
  p=roots(a);
  z=roots(b);
else
  p=a;
  z=b;
end;

hold on;
plot(real(z),imag(z),zc);
plot(real(p),imag(p),pc);
ax=axis;
ax(1)=ax(1)-0.1*abs(ax(1));
ax(3)=ax(3)-0.1*abs(ax(3));
ax(2)=ax(2)+0.1*abs(ax(2));
ax(4)=ax(4)+0.1*abs(ax(4));
axis(ax);
axis('equal');
t=0:(pi/20):(2*pi);
plot(sin(t),cos(t),'b');
hold off;
xlabel('real');
ylabel('imag');
grid;
title(gtitle);

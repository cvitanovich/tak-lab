function e=winenerg(h,wl,n1,n2)
% WINENERG Normalized signal energy inside a rectangular window
%          e=winenerg(h,wl,n1,n2)
%          A rectangular window of length wl is shifted over the
%          impulse resonse h and the energy of the signal inside
%          the window is computed for each position of the window.
%          Optional arguments:
%           n1 --> start index for window (default 1)
%           n2 --> stop index for window (default length(h)-wl)

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 28.09.1994 A. Raab
% (m) 28.07.1995 K. Hartung 




if nargin<4
  n2=length(h)-wl;
  if nargin<3
    n1=1;
  end;
end;

L=n2-n1+1;
e=zeros(L,1);

en=h.^2;   % new K.H.

for n=n1:n2              % compute energy in window
	%  e(n-n1+1)=energy(h(n:n+wl));
	e(n-n1+1)=sum(en(n:n+wl));  % modified K.H.
end;

e=e./max(e);             % normalize energy

function hi=lsqinv(h,N,L,wlp)
% LSQINV Least square inverse filter design
%        hi=LSQINV(h,N,L,wlp)
%        input parameters:
%          h      --> measured impulse response
%          N      --> number of recognized samples (default: length of h -1)
%          L      --> time shift for desired output (default 0)
%          wlp    --> wlp=1 (default) performs a
%                     lowpass filtering of the output signal
%        output parameters:
%          hi --> impulse response of inverse filter
%
% LSQINV solves the set of linear equations
%
%  |rhh(0) rhh(1)   ... rhh(N)  ||hi(0)|   |h(L)  |
%  |rhh(1) rhh(0)   ... rhh(N-1)||hi(1)|   |h(L-1)|
%  | .      .            .      || .   | = | .    |
%  | .      .            .      || .   |   | .    |
%  |rhh(N) rhh(N-1) ... rhh(0)  ||hi(N)|   |h(L-N)|
%
% where h(k)=0 for k<0 and rhh(k) ist the autocorrelation function of h

% References:   Proakis & Manolakis
%               Digital Signal Processing
%               2nd ed.,
%               Macmillan Publishing Company,
%               New York, 1992, p.663
%
%              A. Raab
%              Untersuchungen zum Entwurf von Filtern zur
%              Approximation von Aussenohruebertragungsfunktionen
%              Diploma thesis, Lehrstuhl fuer allgemeine Elektrotechnik
%              und Akustik, Ruhr Universitaet Bochum, 1994

if nargin<4
  wlp=1;              % do lowpass filtering
  if nargin<3
    L=0;              % time offset delta(n-L)
    if nargin<2
      N=length(h)-1;  % inverse filter length
    end;
  end;
end;

L=L+1;                % shift index to keep MATLAB notation
N=N+1;

h=h(1:N); h=h(:);
rhh=xcorr(h);         % compute autocorrelation function
rhh=rhh(N:length(rhh));
R=toeplitz(rhh);      % set up right hand side of equation

hv=flipud([zeros(N-L,1);h(1:L)]); % left hand side of equation
hinv=R\hv;            % solve equation

if (wlp==1)
  % load lp15k;
  b=fir1(41,(16/22.05));
  hi=conv(b(:),hinv(:));    % lowpass filtering
  % hi=filter(b,1,hinv);
 else
  hi=hinv;
end;

hi=hi(:);

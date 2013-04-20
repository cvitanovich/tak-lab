function b=lsqinv2(h,L,Wn,Fs)
% LSQINV2 Least square inverse filter design with low pass or band pass filtering
%        b=LSQINV2(h,L,Wn,Fs)
%        input parameters:
%          h      --> measured impulse response
%          L      --> length of desired impulse response
%	   Wn     --> cutoff frequency in Hz of low-pass filter 0< Wn < Nyquist/1.1
%                 default: no filtering
%	   Fs     --> sampling frequency in Hz (Default Fs = 2)
%        output parameters:
%          b --> impulse response of inverse filter
%
% LSQINV2 solves the set of linear equations
%
%  |rhh(0) rhh(1)   ... rhh(N)  ||b(0)|   |rdh(0)  |
%  |rhh(1) rhh(0)   ... rhh(N-1)||b(1)|   |rdh(1)|
%  | .      .            .      || .   | = | .    |
%  | .      .            .      || .   |   | .    |
%  |rhh(N) rhh(N-1) ... rhh(0)  ||b(N)|   |rdh(N)|
%
% where h(k)=0 for k<0 and rhh(k) is the autocorrelation function of h

% References:   Proakis & Manolakis
%               Digital Signal Processing
%               2nd ed.
%
%              A. Raab  modified by Klaus Hartung
%
% low pass filtering response bandpass filtering
% before linear equations are solved
% solving rdh(n) instead of h(1)


wlp=1;
if (nargin<3)
  wlp=0;              % no lowpass filtering
else
  wlp=1;
  rp=0.1;
  rs=150;
%***************** 
  f=[Wn 1.1*Wn];
%*****************
  m=[1 0];
  dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
  [n,fo,mo,w]=remezord(f,m,dev,Fs);
  d_=remez(n+30,fo,mo,w);
end;

if (nargin < 2)   L=length(h); end;    % time offset delta(n-L)

if (L < n) error('more coeff'); end;

d= zeros(1,L);

h=h(:);
h_ = zeros(1,L);
h_(1:length(h))=h;
h=h_(:);
rhh=xcorr(h);
rhh=rhh(L:length(rhh));

if (wlp ==1)
  d(1:length(d_))=d_;
else
  d=zeros(1,L);
  d(1)=1;
end;

d=d(:);
[m,n]=size(d);

rdh=xcorr(d,flipud(h));			%cross correlation function
rdh=rdh(m:length(rdh));
R=toeplitz(rhh);			% set up right side of eqn
size(R)
size(rdh)
b=R\rdh;				% solve eqn

return
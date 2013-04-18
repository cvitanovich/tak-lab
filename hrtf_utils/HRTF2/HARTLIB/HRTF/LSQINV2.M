function b=lsqinv2(h,L,Wn, Fs)
% LSQINV2 Least square inverse filter design with lowpass or bandpass-filtering 
%        hi=LSQINV2(h,L,Wn,Fs)
%        input parameters:
%          h      --> measured impulse response
%          L      --> length of desired impulse response
%          Wn     --> cutoff frequency of low-pass filter in Hz 
%                     default: no filtering
%	   Fs     --> sampling rate in Hz (default Fs = 2)
%        output parameters:
%          b --> impulse response of inverse filter
%
% LSQINV solves the set of linear equations
%
%  |rhh(0) rhh(1)   ... rhh(N)  ||b(0)|   |rdh(0)|
%  |rhh(1) rhh(0)   ... rhh(N-1)||b(1)|   |rdh(1)|
%  | .      .            .      || .  | = | .    |
%  | .      .            .      || .  |   | .    |
%  |rhh(N) rhh(N-1) ... rhh(0)  ||b(N)|   |rdh(N)|
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


% modified by Klaus Hartung
% Lehrstuhl fuer allgemeine Elektrotechnik
%              und Akustik, Ruhr Universitaet Bochum, 1995
% low pass filtering resp. bandpass filtering before linear equations are solved.
% solving rdh(n) instead of h(l)


if (nargin < 2) error('Not enough arguments'); end;

wlp=1;



if (nargin<3)
  wlp=0;              % no filtering
else
	if (nargin < 4) Fs=2; end;
	wlp = 1;
	%*****************************
	rp =  0.1;  % dB Ripples in passband
  	rs = 60;    % dB Attenuation

	   f=[Wls Wlp Wup Wus]; % corner frequencies
        m_dB=[ -60   0   -60   ]; % dB magnitude in bands
      dev_dB=[   3   0.1   3   ]; % DB deviation in bands
	%*****************************
 	 m = 10.^(m_dB/20);
         dev =10.^(m_dev/20);
%	  dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
  	[n,fo,mo,w]=remezord(f,m,dev,Fs);
 	d_=remez(n,fo,mo,w);
end;

if (nargin<2)
    L=length(h);              % time offset delta(n-L)
end;

if (L < n) error('more coeff'); end; 

d=zeros(1,L);

h=h(:);
h_ = zeros(1,L);
h_(1:length(h))=h;
h=h_(:);
rhh=xcorr(h);                     % compute autocorrelation function
rhh=rhh(L:length(rhh));

if (wlp==1) 
  d(1:length(d_))=d_;
else
  d=zeros(1,L);
  d(1)=1;
end;


d=d(:);
[m,n]=size(d);

rdh=xcorr(d,flipud(h));           % cross-correlation function
rdh=rdh(m:length(rdh));
R=toeplitz(rhh);                  % set up right hand side of equation
size(R)
size(rdh)
b=R\rdh;                          % solve equation

return


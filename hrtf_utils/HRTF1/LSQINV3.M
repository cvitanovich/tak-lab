function b=lsqinv3(h,L,Wn,Fs,gf,gain)
% LSQINV3 Least square inverse filter design with low pass or band pass filtering
%        b=LSQINV3(h,L,Wn,Fs,gf, gain)
%        input parameters:
%          h      --> measured impulse response
%          L      --> length of desired impulse response
%	   Wn     --> corner frequencies of bp filter in Hz
%		[lower_stop_band, lower_pass_band, upper_pass_band, upper_stop_band]
%
%			default: no filtering
%	   Fs     --> sampling frequency in Hz (Default Fs = 30000)
%	   gf     --> plot of transfer functions and impulse responses
%		(1: yes, 2: no, default no)
%
%     gain dB attenuation at each band of Wn:
%		[lower_stop_band, lower_pass_band, passband,passband, upper_pass_band, upper_stop_band]
%
%        output parameters:
%          b --> impulse response of inverse filter
%
% LSQINV solves the set of linear equations
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

[L_h,n_channel] = size(h);
if (nargin < 1) error('ERROR IN LSQINV3: no impulse response specified'); end;
if (nargin < 2) L=L_h; end;
if (nargin < 4) Fs=30000; end;
if (nargin < 5) gf=0; end;
if (nargin < 6)
     gain = [-40  -40  0  0  -60  -60];                 % dB mag in bands
%  gain = [-100  -100  0  0  -60  -60];                 % dB mag in bands
end

if (L < L_h) error('ERROR IN LSQINV3: length too short');end;

d = zeros(1,L);
h = [h; zeros(L-L_h,n_channel)];
b = zeros(L,n_channel);


if (nargin<3)
  d(1) = 1;              % no lowpass filtering
else

  f = [ 0 Wn(1) Wn(2) Wn(3) Wn(4) Fs/2];             % corner freqs
  m_dB = gain;
  m = 10.^(m_dB/20);
  d = firls(L-1, 2*f/Fs, m);
  d = d(:);
  if (gf==1)
    D=rfft(d);
    irplot(d,Fs,'Target Band Pass IR');grcntrl;
    magplot(D,Fs, 'Target Band Pass TF - Log Mag','',Wn(1),Wn(4));grcntrl;
    phasplot(D,Fs, 'Target Band Pass TF - Phase ','',Wn(1),Wn(4));grcntrl;
    tgrpplot(D,Fs,'Target Band Pass TF - Group Delay','',Wn(1),Wn(4));grcntrl;
  end;
end;

for nc=1:n_channel
  rhh=xcorr(h(:,nc));
  rhh=rhh(L:length(rhh));
  [m,n]=size(d);
  rdh=xcorr(d,flipud(h(:,nc)));
  rdh=rdh(m:length(rdh));
  R=toeplitz(rhh);
  b(:,nc)=R\rdh;
 
 if (gf==1)
    irplot(b,Fs,'Inverse System IR');grcntrl;
    B=rfft(b);
    magplot(B,Fs, 'Inverse System TF - Log Mag','',Wn(1),Wn(4));grcntrl;
    phasplot(B,Fs, 'Inverse System TF - Phase ','',Wn(1),Wn(4));grcntrl;
    tgrpplot(B,Fs,'Inverse System TF - Group Delay','',Wn(1),Wn(4));grcntrl;
  end;
end;
return
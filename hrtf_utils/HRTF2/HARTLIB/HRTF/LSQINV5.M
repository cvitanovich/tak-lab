function b=lsqinv5(h,L,Wn, Fs,gf)
% LSQINV5 Least square inverse filter design with lowpass or bandpass-filtering 
%
%        hi=LSQINV3(h,L,Wn,Fs)
%        input parameters:
%          h      --> measured impulse response
%          L      --> length of desired impulse response
%          Wn     --> corner frequencies of bp filter in Hz 
%		     [lower_stop_band, lower_pass_band, upper_pass_band, upper_stop_band] 
%                     default: no filtering
%	   Fs     --> sampling rate in Hz (default Fs = 2)
%	   gf     --> plot of  transfer functions and impulse responses 
%		      (1: yes 0: no, default no);
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

[L_h,n_channel]=size(h);
if (nargin<1) error('no impulse response specified'); end; 
if (nargin<2)  L=L_h;  end;            
if (nargin<4)  Fs=44100; end;
if (nargin<5)  gf=0; end;

d=zeros(1,L); 
z=zeros(L-L_h,n_channel);
h=[h;z];
b=zeros(L,n_channel);


if (nargin<3)  
  d(1)=1;
 else	
	Wls=Wn(1);
	Wlp=Wn(2);	
	Wup=Wn(3);
        Wus=Wn(4);
%	Wls=0.5*Wlp;
%	Wus=1.03*Wup; if (Wup>Fs) Wup=0.99*Fs; end;
Rp=0.5;
Rs=60;
Wp=10500;
Ws=13500;
devs = [ (10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20) ];
[n,Fo,Mo,w]= remezord([Wp Ws],[1 0],devs,Fs);
n=L-1;
d=remez(n,Fo,Mo,w);
size(d)
d=[d  zeros(1,L-length(d))];
size(d)
pause
%	f    = [  0   Wls  Wlp  Wup  Wus   Fs/2 ]; % corner frequencies
%	m_dB = [-40   -40   0    0   -60  -60 ]; % dB magnitude in bands
%	m = 10.^(m_dB/20);
%	d=firls(L-1, 2*f/Fs, m);
	d=d(:);
	if (gf==1) 
		D=rfft(d);
		irplot(d,Fs,'Target Band Pass IR');grcntrl;
		magplot(D,Fs, 'Target Band Pass TF - Log Mag','',Wls,Wus);grcntrl;
		%phasplot(D,Fs, 'Target Band Pass TF - Phase ','',Wls,Wus);grcntrl;
		%tgrpplot(D,Fs, 'Target Band Pass TF - Group Delay','',Wls,Wus);grcntrl;
pause
	end;
end;

for nc=1:n_channel

	rhh=xcorr(h(:,nc));                     % compute autocorrelation function
	rhh=rhh(L:length(rhh));



	[m,n]=size(d);
	rdh=xcorr(d,flipud(h(:,nc)));           % cross-correlation function
	rdh=rdh(m:length(rdh));
	R=toeplitz(rhh);                  % set up right hand side of equation
	b(:,nc)=R\rdh;                          % solve equation

	if (gf==1)
		irplot(b,Fs, 'Inverse System IR');grcntrl;
		B=rfft(b);
		magplot(B,Fs, 'Inverse System TF - Log Mag','',Wls,Wus), grcntrl;
		phasplot(B,Fs, 'Inverse System TF - Phase','',Wls,Wus),grcntrl;
		tgrpplot(B,Fs, 'Inverse System TF - Group Delay','',Wls,Wus), grcntrl;
	end

end
return


% end of file lsqinv3.m
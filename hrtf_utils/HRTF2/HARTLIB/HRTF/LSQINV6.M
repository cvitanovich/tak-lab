function b=lsqinv6(h, Wn, Fs,gf)
% LSQINV6 Least square inverse filter design with lowpass or bandpass-filtering 
%
%        hi=LSQINV6(h, Wn,Fs, gf)
%        input parameters:
%          h      --> measured impulse response
%%         L      --> length of desired impulse response
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

%d=zeros(1,L); 
%z=zeros(L-L_h,n_channel);
%h=[h;z];
%b=zeros(L,n_channel);


d=zeros(L_h,1);
d(1)=1;
d=shift(d,floor(L_h/2));

if (nargin<2)
   d2=zeros(L,1);  
   d2(1)=1;
   [m,p]=max(d);
   L_d=length(d);
   m_shift=floor(L_d/2)-m
%   d2=shift(d2,m_shift);
 else	
   Rsl=10;
   Rp=0.5;
   Rsh=80;
  Wp = Wn(1);
  Ws = Wn(2);
%   Wp=13000;
%   Ws=15500;
   devs = [ (10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rsh/20) ];
%   devs = [ 10^(-Rsl/20) (10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rsh/20) ];

%               
   [n,Fo,Mo,w]= remezord([Wp Ws],[1 0],devs,Fs);
%    [n,Fo,Mo,w]=remezord([ 300 500 Wp Ws],[.001,1,0],devs,Fs);
%n=32
   d2=remez(n,Fo,Mo,w);


%    d=remez(n,[0 Wp Ws 22050]/22050,[1 1 1e-5 1e-5]);
%    d_mean=mean(d)
%    d=d-d_mean*zeros(size(d));
%    d=[d zeros(length(L-d,1))];     
%   size(d)

   d=d(:);
	if (gf==1) 
		D=rfft(d);
%		irplot(d,Fs,'Target Band Pass IR');grcntrl;
%		magplot(D,Fs, 'Target Band Pass TF - Log Mag','',100,Fs/2);grcntrl;
		%phasplot(D,Fs, 'Target Band Pass TF - Phase ','',Wls,Wus);grcntrl;
		%tgrpplot(D,Fs, 'Target Band Pass TF - Group Delay','',Wls,Wus);grcntrl;
                pause
	end;
end; % of if nargin <2


d=lsqfirin(h,d,0);

for k=1:n_channel
  b(:,k)=conv(d(:,k),d2);
end
return


% end of file lsqinv3.m

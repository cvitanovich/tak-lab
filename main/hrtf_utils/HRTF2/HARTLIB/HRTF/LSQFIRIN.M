function b=lsqfirin(h,d,gf)
% LSQFIRIN Least square inverse filter design with lowpass or bandpass-filtering 
%
%        b=LSQFIRIN(h,d)
%        input parameters:
%          h      --> measured impulse response
%          d     --> desired output sequence
%         gf     --> grafic on / off		
%		      (1: yes 0: no, default no);

%        output parameters:
%          b --> impulse response of inverse filter
%
% LSQFIRIN solves the set of linear equations
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

% rewritten based on lsqinv3... 20.7.96 Klaus Hartung

% solving rdh(n) instead of h(l)

Fs=1;
Wls=0;
Wus=Fs/2;
if (nargin<1) error('no impulse response specified'); end; 
if (nargin<2) error('no desired output sequence specified'); end;
if (nargin<3) gf=0; end;
           
[L_h,n_channel]=size(h);
d=d(:);
L_d=length(d);
%if (L_h > L_d) error('Desired output sequence too short'); end;


L_max=max([L_h,L_d]);
L_min=min([L_h,L_d]);

if (L_h~=L_d) 
  L_max=max([L_d,L_h]);
end;

if (L_h<L_max)
  h=[h; zeros(L_max-L_h,1)];
end

%h_scale=1/max(h);
%d_scale=1/max(d);
%d=d*d_scale;
%h=h*h_scale;
%irplot(h,44100);
%irplot(d,44100);
for nc=1:n_channel
disp(['Computing channel ' int2str(nc) ' of ' int2str(n_channel)]);
	rhh=xcorr(h(:,nc));                     % compute autocorrelation function
	rhh=rhh(L_max:length(rhh));
	[m,n]=size(d);
        rdh=xcorr(d,h(:,nc));
        lrdh=length(rdh);
        lrdh2=ceil(lrdh/2);
	rdh=flipud(rdh);
	rdh=rdh(lrdh2:lrdh);
%        rdh=rdh/max(rdh);
%      	 rhh=rhh/max(rhh);
	R=toeplitz(rhh);                        % set up right hand side of equation
	b(:,nc)=R\rdh;                          % solve equation

%       b(lrdh2/2:lrdh2,nc)=zeros(lrdh2/2+1,1);
	if (gf==1)
		irplot(b,Fs, 'Inverse System IR');grcntrl;
		%B=rfft(b);
                u=conv(b(:,nc),h(:,nc));
                U=rfft(u,4096*2);
                Wls=100;
                Wus=20000;
		irplot(u,44100, 'Desired Sysem Inpulse response');

 dd=[d; zeros(length(u)-length(d),1)];	
 irplot((u-dd),44100, 'Error');
		magplot(U,44100, 'Desired TF - Log Mag','',Wls,Wus);
%		phasplot(U,44100, 'Inverse System TF - Phase','',Wls,Wus),grcntrl;
%		tgrpplot(U,44100, 'Inverse System TF - Group Delay','',Wls,Wus), grcntrl;
	end
%clc

end                             % of foor loop



% end of file lsqinv3.m
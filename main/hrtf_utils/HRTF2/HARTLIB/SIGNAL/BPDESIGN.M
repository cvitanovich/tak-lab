function bpdesign
Wls=100;
Wlp=200;
Wup=16000;
Wus=17000;
dev=[];
Fs = 44100;

%	rp =  1;  % dB Ripples in passband
%  	rs = 40;    % dB Attenuation
%
	   f=[0     Wls  Wlp  Wup  Wus   Fs ]; % corner frequencies
	   m=[0      0    1    1     0    0];
        m_dB=[-40   -40   0    0    -50  -50 ]; % dB magnitude in bands
%      dev_dB=[ -80         3          -60  ]; % DB deviation in bands
	%*****************************
%dev=[10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]
 	 m = 10.^(m_dB/20)
%        dev =(10.^(dev_dB/20))
%          dev = [0.0001  0.1  0.0001]
%	  dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
%  	[n,fo,mo,w]=remezord(f,m,dev,Fs)
% 	d_=remez(n,fo,mo,w);
d_=firls(512,f/Fs,m);
D=rfft(d_,4096);
magplot(D,Fs);
grcntrl;
irplot(d_',Fs);
grcntrl

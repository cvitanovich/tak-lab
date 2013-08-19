function [p,delay]=value2place(itd, n_taps, center_delay)

% VALUE2PLACE codes scalar values in vector pattern
%
% [p,delay]=VALUE2PLACE(v, n, c)
%
% v: scaler value
% n: number of taps
% c: center position on vector

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 16.09.1996

% Calls         Toolbox
radbas		nn/nn

[m,n]=size(itd);

b=n_taps;

delta=2/(n_taps-1);

e=max(max(abs(itd)))

min_delay=-e;
max_delay=e;

for f=1:m
  delay(f,:)=linspace(min_delay(f),max_delay(f),n_taps);
end

p=[];
pd=[];

for d=1:n
	for f=1:m
		t = itd(f,d)*ones(1,n_taps);
                dt=t-delay(f,:);
		excitation=radbas(dt,b);
		pd=[pd; excitation];
	end;

	[x,y]=size(pd);
        delay;
	pd=reshape(pd,x*y,1);
	p=[p pd];


end

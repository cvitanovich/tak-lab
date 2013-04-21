function fcoefs=make1_ERBFiltA(Fs,cF)
% function [fcoefs]=make1_ERBFilt1(Fs,cF);
% This function computes the filter coefficients for a single 
% Gammatone filter.  These filters were defined by Patterson and 
% Holdworth for simulating the cochlea.  
% 
% The result is returned as an array of filter coefficients.  Each row 
% of the filter arrays contains the coefficients for four second order 
% filters.  The transfer function for these four filters share the same
% denominator (poles) but have different numerators (zeros).  All of these
% coefficients are assembled into one vector that the ERBFilterBank 
% can take apart to implement the filter.
%
%argins:
% Fs is the sampling frequency
% cF is the center frequency
%
% NOTE:! this implementation fixes a problem in the original code by
% computing four separate second order filters.  This avoids a big
% problem with round off errors in cases of very small cfs (100Hz) and
% large sample rates (44kHz).  The problem is caused by roundoff error
% when a number of poles are combined, all very close to the unit
% circle.  Small errors in the eigth order coefficient, are multiplied
% when the eigth root is taken to give the pole location.  These small
% errors lead to poles outside the unit circle and instability.  Thanks
% to Julius Smith for leading me to the proper explanation.

% Execute the following code to evaluate the frequency
% response of a 10 channel filterbank.
%	fcoefs = MakeERBFilters(16000,10,100);
%	y = ERBFilterBank([1 zeros(1,511)], fcoefs);
%	resp = 20*log10(abs(fft(y')));
%	freqScale = (0:511)/512*16000;
%	semilogx(freqScale(1:255),resp(1:255,:));
%	axis([100 16000 -60 0])
%	xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');

% Rewritten by Malcolm Slaney@Interval.  June 11, 1998.
% (c) 1998 Interval Research Corporation  

% Version to implement owl's bandwidths...

if nargin <2
   error('not enough input arguments')
   return
end

% Change the following three parameters if you wish to use a different
% ERB scale.  Must change in ERBSpace too.
%EarQ = 9.26449;				%  Glasberg and Moore Parameters
%minBW = 24.7;
%order = 1;

%ERB = ((cF/EarQ).^order + minBW^order).^(1/order);


T=1/Fs;
EarQ = cF / (2*50); %Make narrower gammatone filter -- NOT the owl's auditory nerve!
%EarQ = 1.9 * 0.074 * (cF .^ .504);
ERB = cF ./ EarQ;
B=1.019*2*pi*ERB;

A0 = T;
A2 = 0;
B0 = 1;
B1 = -2*cos(2*cF*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*cF*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cF*pi*T)./ ...
		exp(B*T))/2;
A12 = -(2*T*cos(2*cF*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cF*pi*T)./ ...
		exp(B*T))/2;
A13 = -(2*T*cos(2*cF*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cF*pi*T)./ ...
		exp(B*T))/2;
A14 = -(2*T*cos(2*cF*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cF*pi*T)./ ...
		exp(B*T))/2;

gain = abs((-2*exp(4*i*cF*pi*T)*T + ...
                 2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
                         (cos(2*cF*pi*T) - sqrt(3 - 2^(3/2))* ...
                          sin(2*cF*pi*T))) .* ...
           (-2*exp(4*i*cF*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
              (cos(2*cF*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cF*pi*T))).* ...
           (-2*exp(4*i*cF*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
              (cos(2*cF*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cF*pi*T))) .* ...
           (-2*exp(4*i*cF*pi*T)*T + 2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
           (cos(2*cF*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cF*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4*i*cF*pi*T) +  ...
           2*(1 + exp(4*i*cF*pi*T))./exp(B*T)).^4);
	
allfilts = ones(length(cF),1);
fcoefs = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];

if (0)						% Test Code
	A0  = fcoefs(:,1);
	A11 = fcoefs(:,2);
	A12 = fcoefs(:,3);
	A13 = fcoefs(:,4);
	A14 = fcoefs(:,5);
	A2  = fcoefs(:,6);
	B0  = fcoefs(:,7);
	B1  = fcoefs(:,8);
	B2  = fcoefs(:,9);
	gain= fcoefs(:,10);	
	chan=1;
	x = [1 zeros(1, 511)];
	y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
		A2(chan)/gain(chan)],[B0(chan) B1(chan) B2(chan)], x);
	y2=filter([A0(chan) A12(chan) A2(chan)], ...
			[B0(chan) B1(chan) B2(chan)], y1);
	y3=filter([A0(chan) A13(chan) A2(chan)], ...
			[B0(chan) B1(chan) B2(chan)], y2);
	y4=filter([A0(chan) A14(chan) A2(chan)], ...
			[B0(chan) B1(chan) B2(chan)], y3);
	semilogx((0:(length(x)-1))*(Fs/length(x)),20*log10(abs(fft(y4))));
end

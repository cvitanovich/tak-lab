function [forward, feedback, cf]=ERBfilt2(fs,cf)
% [forward, feedback, cf]=MakeERBFilters(fs,numChannels) computes the
% filter coefficients for a bank of Gammatone filters.  These
% filters were defined by Patterson and Holdworth for simulating 
% the cochlea.  The results are returned as arrays of filter
% coefficients.  Each row of the filter arrays (forward and feedback)
% can be passed to the MatLab "filter" function, or you can do all
% the filtering at once with the filtBank() function.
%
% The filter bank contains "numChannels" channels that extend from
% half the sampling rate (fs) to "lowFreq".
% if numChannels == 0 then 1/6 octave spacing is used 
% (20 channels);
% otherwise Kopple et al.'s 1993 spacing is used
%
%
% MODIFIED TO DO ONLY ONE FILTER at cf

T=1/fs;
                         

% 2.5 * reflects change from Q10dB 
% to something near critical band (ca. Q4dB ??)

EarQ = 2.5 * 0.074 * (cf .^ .504);
ERB = cf ./ EarQ;


B=1.019*2*pi*ERB;
gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
				2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
				(cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
				sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cf*pi*T))).* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
           (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
           2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);
feedback=zeros(length(cf),9);
forward=zeros(length(cf),5);
forward(:,1) = T^4 ./ gain;
forward(:,2) = -4*T^4*cos(2*cf*pi*T)./exp(B*T)./gain;
forward(:,3) = 6*T^4*cos(4*cf*pi*T)./exp(2*B*T)./gain;
forward(:,4) = -4*T^4*cos(6*cf*pi*T)./exp(3*B*T)./gain;
forward(:,5) = T^4*cos(8*cf*pi*T)./exp(4*B*T)./gain;
feedback(:,1) = ones(length(cf),1);
feedback(:,2) = -8*cos(2*cf*pi*T)./exp(B*T);
feedback(:,3) = 4*(4 + 3*cos(4*cf*pi*T))./exp(2*B*T);
feedback(:,4) = -8*(6*cos(2*cf*pi*T) + cos(6*cf*pi*T))./exp(3*B*T);
feedback(:,5) = 2*(18 + 16*cos(4*cf*pi*T) + cos(8*cf*pi*T))./exp(4*B*T);
feedback(:,6) = -8*(6*cos(2*cf*pi*T) + cos(6*cf*pi*T))./exp(5*B*T);
feedback(:,7) = 4*(4 + 3*cos(4*cf*pi*T))./exp(6*B*T);
feedback(:,8) = -8*cos(2*cf*pi*T)./exp(7*B*T);
feedback(:,9) = exp(-8*B*T);
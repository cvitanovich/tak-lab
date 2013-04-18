function hz = bark2hz(bark)

% function hz = bark2hz(bark)
%
% BARK2HZ changes frequency value from the nonlinear 'bark' - scale to the 
%         linear frequency scale in 'Hz'.
%         The formular is taken from: Sottek R. (1993); Modelle zur Signal-
%         verarbeitung im menschlichen Gehör; Dissertation am Lehrstuhl fuer
%         Elektrotechnik RWTH Aachen.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 23.09.1994 Joerg Schmatter
% rewritten by Klaus Hartung 26.7.95



r=find(bark>24);
if (length(r) ~= 0) 
 disp('bark values are greater than 24')
 hz=-1;
 return;
end;

r=find(bark<0);
if (length(r) ~= 0) 
 disp('bark values are less than 0')
 return;
 hz=-1;

end;

nom = bark + 0.53;
denom=26.28 - bark;

hz=(nom./denom).*1960;

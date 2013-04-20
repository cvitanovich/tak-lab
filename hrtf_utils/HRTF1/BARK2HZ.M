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

hz = 1960 * ( (bark + 0.53)/(26.28 - bark) );

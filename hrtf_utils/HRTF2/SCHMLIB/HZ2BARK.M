function bark = hz2bark(hz)

% function bark = hz2bark(hz)
% 
% BARK changes frequency value from the linear 'Hz' - scale to the nonlinear
%      frequency scale in 'bark'.
%      The formular is taken from: Sottek R. (1993); Modelle zur Signal-
%      verarbeitung im menschlichen Gehör; Dissertation am Lehrstuhl fuer
%      Elektrotechnik RWTH Aachen.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


bark = 13.0 * atan(0.76/1000.0*hz) + 3.5 * atan((hz/7500.0)^2);
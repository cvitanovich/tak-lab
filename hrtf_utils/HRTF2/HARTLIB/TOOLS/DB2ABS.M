function a=db2abs(dbv)


% DB2ABS converts value in dB to abs value
%
%function a=db2abs(dbv)

% Klaus Hartung
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 12.5.1996

a=10.^(dbv/20);

return


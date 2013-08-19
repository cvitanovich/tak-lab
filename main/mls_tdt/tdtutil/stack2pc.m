function [var_name] = stack2pc()
%[var_name] = stack2pc(Fs)
%var_name:		must name a time-domain real signal
%Fs:			sampling frequency
%


while (S232('APactive'))  pause(0); end
var_name = S232('pop16');

return
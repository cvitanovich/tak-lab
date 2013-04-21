function [damabuf] = pc2dama(var_name,Fs,scalefact,DAMA_number)
%pc2dama(var_name,Fs,scalefact,DAMA_number):		Go from the variable space into an int16 DAMA buffer
%var_name:		must name a time-domain real signal
%Fs:			sampling frequency
%scalefact:		specify scaling
%

if(~isreal(var_name))
   error('Sound variable must contain a real array');
else
   [r,c] = size(var_name);
   if(r > 1)
      var_name = var_name';
      disp('WARNING: Variable must be a row vector');
   end
end

%Set mean amplitude to 0
mean_var_name = mean(var_name);
var_name = var_name - mean_var_name;

%Apply scaling
if(nargin < 3 & (max(abs(var_name)) > 0))
   var_name = round(10000 * var_name/max(abs(var_name)));
elseif(nargin >= 3)
   var_name = round(scalefact * var_name);
end

%Set mean amplitute to 0
mean_var_name = mean(var_name);
var_name = var_name - mean_var_name;

%Make the DAMA buffer
while (S232('APactive') == 1)  pause(0); end
S232('push16',var_name,length(var_name));       % push onto stack

if(nargin < 4)
   while (S232('APactive') == 1)  pause(0); end
   DAMA_number = S232('_allot16',length(var_name));
end

while (S232('APactive') == 1)  pause(0); end
S232('qpop16',DAMA_number);                     % pop from stack to DAMA

return
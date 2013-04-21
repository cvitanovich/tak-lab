function [] = pc2stack(var_name)
%pc2stack(var_name,Fs):		Go from the variable space into a DAMA buffer
%var_name:		must name a time-domain real signal
%Fs:			sampling frequency
%

if(~isreal(var_name))
   error('Variable must contain a real array');
else
   [r,c] = size(var_name);
   if(r > 1)
      var_name = var_name';
      disp('WARNING: Variable must be a row vector');
   end
end

mean_var_name = mean(var_name);
var_name = var_name - mean_var_name;
var_name = round(10000 * var_name/max(abs(var_name)));
mean_var_name = mean(var_name);
var_name = round(var_name - mean_var_name);

while (S232('APactive'))  pause(0); end
S232('push16',var_name,length(var_name));

return
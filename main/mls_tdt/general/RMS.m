function [rmsval] = RMS(var_name)
%[rmsval] = RMS(var_name)

rmsval = sqrt(mean(var_name.^2));

return;
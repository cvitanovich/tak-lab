function y = filtbank(forward,feedback,x);

% y = ERBFilterBank(forward,feedback,x)
% This function filters the waveform x with the array of filters
% specified by the forward and feedback paarameters.  Each row
% of the forward and feedback parameters are the parameters 
% to the MATLAB builtin function "filter"

[rowf, colf] = size(feedback);
[rowx,colx] = size(x);
if (rowx > colx) 
   x = x';
   [rowx,colx] = size(x);
end
   
y = zeros(rowf, colx);
for i = 1:rowf
   y(i,:) = filter(forward(i,:),feedback(i,:),x);
end

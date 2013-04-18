function [u_s_signal,new_xaxis,newlines] = u_sample(signal,...
                               xaxis,n_directions,n_lines,resample_factor)

% function [u_s_signal,new_xaxis,newlines] = u_sample(signal,...
%                               xaxis,n_directions,n_lines,resample_factor)
%
% U_SAMPLE resamples the in the matrix 'signal' included datas with a lower
%          sampling rate after lowpass filtering. The datas are resampled at
%          1/resample_factor times of the original sample rate. The resulting
%          resampled matrix 'u_s_signal' is 'resample_factor' times shorter.
%          The number of lines of the resampled datas is given in 'newlines'
%          and the vector including the frequencies of sampling is adapted
%          to the new signal matrix. 
%
%          For further information see: decimate
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

u_s_signal = [];

for i=1:n_directions
  sig = decimate(signal(i,:),resample_factor);
  u_s_signal = [u_s_signal;sig];
end;

[rows,newlines] = size(u_s_signal);
new_xaxis = [newlines];
index = 1;
for i=1:newlines
  new_xaxis(i) = xaxis(index);
  index = index + resample_factor;
end;
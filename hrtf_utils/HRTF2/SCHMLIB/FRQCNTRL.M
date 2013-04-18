function selected_frq = frqcntrl(frq,n_lines,begin_of_signal,end_of_signal,...
                                 xaxis)

% function selected_frq = frqcntrl(frq,n_lines,begin_of_signal,end_of_signal,...
%                                  xaxis)
%
% FRQCNTRL returns the value of the sampled frequency, that is the nearest
%          to the choosen frequency, if the choosen one is not sampled exactly.
%
%          meaning of the input arguments:
%          - frq: choosen frequency
%          - begin_of_signal
%          - end_of_signal
%          - xaxis: vector including the points of sampling
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


% evaluating the number of line for the selected_frq

found = 0;
selected_line = 1;
while (found~=1)
  left_frq = xaxis(selected_line);
  right_frq = xaxis(selected_line+1);
  if (left_frq<=frq)&(right_frq>=frq)
    found = 1;
  else
    selected_line = selected_line + 1;
  end;
end;

% evaluate nearest frequency of the selected frequency

left_delta = frq - left_frq;
right_delta = right_frq - frq;
if (left_delta<=right_delta)
  selected_frq = xaxis(selected_line);
else
  selected_line = selected_line + 1;
  selected_frq = xaxis(selected_line);
end;




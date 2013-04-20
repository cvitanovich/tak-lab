function [selected_frq,selected_line] = ia_s_frq(frq_axis,n_iaural)

% function [selected_frq,selected_line] = ia_s_frq(frq_axis,n_iaural)
%
% IA_S_FRQ selects a frequency for plotting interaural signal differences.
%          The user is asked for the selected frequency, which then is compared
%          with the sampled frequencies. If necessary the selected frequency
%          is changed to the nearest sampling point. The available frequency
%          and the equivalent number of sampling line are returned.
%
%          Description of the input arguments:
%          - frq_axis: vector including the frequencies in [Hz], for which
%                      the interaural differences are evaluated.
%          - n_iaural: number of interaural analysed frequencies
%
%          Description of the output arguments:
%          - selected_frq: really analysed frequency
%          - selected_line: number of line of the 'selected_frq'
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter

disp(['Begin of signal: ',num2str(frq_axis(1))]);
disp(['End of signal: ',num2str(frq_axis(n_iaural))]);

selected_frq = [];
while isempty(selected_frq)      % selecting frequency
  selected_frq = input('Select frequency (in Hz) to analyse: ');
  if selected_frq<0
     selected_frq = -selected_frq;
  end;
  if selected_frq<frq_axis(1)
    selected_frq = frq_axis(1);
  elseif selected_frq>frq_axis(n_iaural)
    selected_frq = frq_axis(n_iaural);
  end;
end;

% evaluating the number of line for the selected_frq

found = 0;
selected_line = 1;
while (found~=1)
  left_frq = frq_axis(selected_line);
  right_frq = frq_axis(selected_line+1);
  if (left_frq<=selected_frq)&(right_frq>=selected_frq)
    found = 1;
  else
    selected_line = selected_line + 1;
  end;
end;

left_delta = selected_frq - left_frq;
right_delta = right_frq - selected_frq;
if (left_delta<=right_delta)
  selected_frq = frq_axis(selected_line);
else
  selected_line = selected_line + 1;
  selected_frq = frq_axis(selected_line);
end;
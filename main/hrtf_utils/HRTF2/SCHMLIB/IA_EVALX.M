function [level_difference,time_difference,frq_axis] = ...
         ia_evalx(first_channel,second_channel,n_directions,...
                  n_lines,sampling_rate,begin_of_signal,end_of_signal,xaxis)

% function [level_difference,time_difference,frq_axis] = ...
%          ia_evalx(first_channel,second_channel,n_directions,...
%                   n_lines,sampling_rate,begin_of_signal,end_of_signal,xaxis)
%
% IA_EVALX evaluates the interaural level and time differences for all
%          directions and all frequencies in a selected range. The frequency
%          range must be selected in [bark].
%
%          meaning of the input arguments:
%          - first_channel: matrix with the absolute magnitude of the signals
%                           first channel of the signals. Each row represents
%                           one direction, each column one sampling line.
%          - second_channel: matrix with the absolute magnitude of the signals
%                            second channel. Each row represents one direction,
%                            each column one sampling line. If the file includes
%                            only one channel, this matrix is empty.
%          - n_directions: number of directions
%          - n_lines: number of lines for each direction and each channel
%          - sampling_rate
%          - begin_of_signal
%          - end_of_signal
%          - xaxis: vector including the points of sampling
%          meaning of the output arguments:
%          - level_difference: matrix including the interaural level
%                              differences. Each row represents one direction,
%                              each column represents the level differences
%                              for a certain frequency.
%          - time_difference: matrix including the interaural time differences.
%                             Each row represents one direction, each column
%                             represents the level differences for a certain
%                             frequency.
%          - frq_axis: vector including the frequencies in [Hz], for which
%                      the interaural differences are evaluated.
%
%          For further information see: leveldif, timediff, earfilt
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


first_bark = hz2bark(begin_of_signal);
last_bark = hz2bark(end_of_signal);

disp(['Begin of signal: ',num2str(first_bark),' bark']);
disp(['End of signal: ',num2str(last_bark),' bark']);

min_bark = [];
while isempty(min_bark)
  min_bark = abs(input('First frequency in ''bark'' ? '));
  if min_bark < first_bark
   min_bark = first_bark;
  end;
  if min_bark==0
    disp('Choose greater value than ''0'' bark');
    min_bark = [];
  end;
end;

max_bark = [];
while isempty(max_bark)
  max_bark = abs(input('Last frequency in ''bark'' ? '));
  if max_bark<=min_bark
    max_bark = [];
  elseif max_bark>last_bark
    max_bark = last_bark;
  end;
end;

disp('Selected range');
disp(['Minimum frequency: ',num2str(min_bark),' bark']);
disp(['Maximum frequency: ',num2str(max_bark),' bark']);

bark_steps = [];
while isempty(bark_steps)
  bark_steps = abs(input('Frequency steps in ''bark'' ? '));
  bark_range = max_bark - min_bark;
  if rem(bark_range,bark_steps)~=0
    bark_steps = [];
    disp('steps must fit in selected bark range');
  end;
end;

% evaluating interaural differences

resolution = (end_of_signal - begin_of_signal) / (n_lines - 1);
level_difference = [];
time_difference = [];
frq_axis = [];
n_iaural = 0;

low_offset_lines = begin_of_signal / resolution;
high_offset_lines = (sampling_rate/2 - end_of_signal) / resolution;
uncutted_number_of_lines = low_offset_lines + n_lines + high_offset_lines;

% for different frequencies

for bark=min_bark:bark_steps:max_bark
  n_iaural = n_iaural + 1;
  frq = bark2hz(bark);
  selected_frq = frqcntrl(frq,n_lines,begin_of_signal,end_of_signal,xaxis);
  frq_axis(n_iaural) = selected_frq;

  H = earfilt(selected_frq,uncutted_number_of_lines,resolution,sampling_rate);

  % for each direction
  for i=1:n_directions
    left_channel = [];
    right_channel = [];
    if low_offset_lines~=0
      left_channel(1:low_offset_lines) = zeros(1,low_offset_lines);
      right_channel(1:low_offset_lines) = zeros(1,low_offset_lines);
    end;
  
    left_channel(low_offset_lines+1:low_offset_lines+n_lines) = first_channel(i,:);
    right_channel(low_offset_lines+1:low_offset_lines+n_lines) = second_channel(i,:);

    if high_offset_lines~=0
      left_channel(low_offset_lines+n_lines+1:low_offset_lines+n_lines+high_offset_lines) = zeros(1,high_offset_lines);
      right_channel(low_offset_lines+n_lines+1:low_offset_lines+n_lines+high_offset_lines) = zeros(1,high_offset_lines);
    end;

    left_X = left_channel .* H;
    right_X = right_channel .* H;

    % evaluate interaural level difference
    level_difference(i,n_iaural) = leveldif(left_X,right_X);

    % evaluate interaural time difference
    time_difference(i,n_iaural) = timediff(selected_frq,uncutted_number_of_lines,sampling_rate,left_X,right_X);
  end;
  disp([num2str(n_iaural),' th of ',num2str(1+(max_bark-min_bark)/bark_steps), ' frequencies evaluated']);
end;
function [level_difference,time_difference,frq_axis] = ...
         ia_eval1(first_channel,second_channel,n_directions,...
                  n_lines,sampling_rate,begin_of_signal,end_of_signal)

% function [level_difference,time_difference,frq_axis] = ...
%          ia_eval1(first_channel,second_channel,n_directions,...
%                   n_lines,sampling_rate,begin_of_signal,end_of_signal)
%
% IA_EVAL1 evaluates the interaural level and time difference for all given
%          directions and a single frequency.
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
%          meaning of the output arguments:
%          - level_difference: column vector including the interaural
%                              level difference for each direction
%          - time_difference: column vector including the interaural
%                             time difference for each direction
%          - frq_axis: including the frequency in [Hz], for which the
%                      interaural differences are evaluated
%
%          For further information see: leveldif, timediff, earfilt
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


min_frq = 100;
max_frq = 15000;

disp(' ');
disp(['Begin of signal: ',num2str(begin_of_signal),' Hz']);
disp(['End of signal: ',num2str(end_of_signal),' Hz']);

% setting range for analysing

selected_frq = [];
while isempty(selected_frq)
  selected_frq = input('Frequency to analyse ? ');
  if selected_frq<begin_of_signal
    selected_frq = begin_of_signal;
  elseif selected_frq>end_of_signal
    selected_frq = end_of_signal;
  end;
end;

if selected_frq<min_frq
  selected_frq = min_frq;
elseif selected_frq>max_frq
  selected_frq = max_frq;
end;

resolution = (end_of_signal - begin_of_signal) / (n_lines-1);
frq_line = round(1 + (selected_frq - begin_of_signal) / resolution);
selected_frq = begin_of_signal + (frq_line - 1) * resolution;

frq_axis(1) = selected_frq;

level_difference = [];
time_difference = [];
low_offset_lines = begin_of_signal / resolution;
high_offset_lines = (sampling_rate/2 - end_of_signal) / resolution;
uncutted_number_of_lines = low_offset_lines + n_lines + high_offset_lines;

H = earfilt(selected_frq,uncutted_number_of_lines,resolution,sampling_rate);

% for each direction
for j=1:n_directions
  left_channel = [];
  right_channel = [];
  if low_offset_lines~=0
    left_channel(1:low_offset_lines) = zeros(1,low_offset_lines);
    right_channel(1:low_offset_lines) = zeros(1,low_offset_lines);
  end;
  
  left_channel(low_offset_lines+1:low_offset_lines+n_lines) = first_channel(j,:);
  right_channel(low_offset_lines+1:low_offset_lines+n_lines) = second_channel(j,:);

  if high_offset_lines~=0
    left_channel(low_offset_lines+n_lines+1:low_offset_lines+n_lines+high_offset_lines) = zeros(1,high_offset_lines);
    right_channel(low_offset_lines+n_lines+1:low_offset_lines+n_lines+high_offset_lines) = zeros(1,high_offset_lines);
  end;

  left_X = left_channel .* H;
  right_X = right_channel .* H;

  % evaluate interaural level difference
  level_difference(j) = leveldif(left_X,right_X);

  % evaluate interaural time difference
  time_difference(j) = timediff(selected_frq,uncutted_number_of_lines,sampling_rate,...
                                left_X,right_X);
  disp([num2str(j),' th direction evaluated']);
end;

level_difference = level_difference';
time_difference = time_difference';
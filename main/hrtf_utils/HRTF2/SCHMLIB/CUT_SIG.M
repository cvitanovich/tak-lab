function [first_channel,second_channel,first_phase,second_phase,n_lines,...         
begin_of_signal,end_of_signal,xaxis]...   =
cut_sig(filetype,first_channel,second_channel,first_phase,second_phase,...            
n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,...             
end_of_signal,xaxis)


% function [first_channel,second_channel,first_phase,second_phase,n_lines,...
%           begin_of_signal,end_of_signal,xaxis]...
%   = cut_sig(filetype,first_channel,second_channel,first_phase,second_phase,...
%             n_directions,n_channels,n_lines,sampling_rate,begin_of_signal,...
%             end_of_signal,xaxis)
%
% CUT_SIG cuts a signal in the selected range if cutting is selected.
%         At first is asked if a range cutting is desired. If not the orignal
%         data are returned. Else the range has to be selected. The reduced
%         matrices are evaluated and returned. The variables 'n_lines',
%         'begin_of_signal', 'end_of_signal' are set with the new values.
%         The vector 'xaxis' is reduced, too.
%
%         For further information about the input and ouput arguments see:
%         READ
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

cutdecision = [];
while isempty(cutdecision)
  cutdecision = input('Whole range to display [y/n] ? ','s');
  if ( (cutdecision~='y')&(cutdecision~='n') )
    cutdecision = [];
  end;
 end;
if (cutdecision=='n')
  if filetype==2
    disp('Select new range in [sec] !');
  elseif ( (filetype>=3)&(filetype<=6) )
    disp('Select new range in [Hz] !');
  elseif filetype==7
    disp('Select new range in [mel] !');
   elseif filetype==8
    disp('Select new range in [bark] !');
   elseif filetype==9
    disp('Select new range in [ERB] !');
   elseif filetype==10
    disp('Select new range in [log10( f / Hz )] !');
   end;
   new_begin = [];
  while isempty(new_begin)
    new_begin = input('New begin of the signal ? ');
    if new_begin<begin_of_signal
      new_begin = begin_of_signal;
    elseif new_begin>=end_of_signal
      new_begin = [];
      disp('Beginning can''t be greater than ending of the signal');
    end;
  end;
  new_end = [];
  while isempty(new_end)
    new_end = input('New end of the signal ? ');
    if new_end>end_of_signal
      new_end = end_of_signal;
    elseif new_end<=new_begin
      new_end = [];
      disp('Ending can''t be smaller than beginning of the signal');
    end;
  end;
   
% changing time or frequency value to number of line
  if filetype==2
    resolution = 1 / sampling_rate;
  else
    resolution = (end_of_signal - begin_of_signal)/ (n_lines - 1);
  end;
  
% changing begin of the signal
  old_begin = begin_of_signal;
  if (new_begin==begin_of_signal)          % first_line is choosen
    first_line = 1;
    begin_of_signal = new_begin;
  else
    delta = new_begin - begin_of_signal;
    first_line = round(delta / resolution) + 1;
    begin_of_signal = begin_of_signal + (first_line-1) * resolution;
  end;
   
% changing end of the signal
  if (new_end==end_of_signal)
    last_line = n_lines;
    end_of_signal = new_end;
  else
    delta = new_end - old_begin;
    last_line = round(delta / resolution) + 1;
    end_of_signal = old_begin + (last_line-1) * resolution;
  end;
  n_lines = 1 + (end_of_signal - begin_of_signal) / resolution;
  disp('selected range: ');
  if filetype==2
    disp(['begin: ',num2str(begin_of_signal),' sec']);
    disp(['end: ',num2str(end_of_signal),' sec']);
  elseif ( (filetype>=3)&(filetype<=6) )
    disp(['begin: ',num2str(begin_of_signal),' Hz']);
    disp(['end: ',num2str(end_of_signal),' Hz']);
  elseif filetype==7
    disp(['begin: ',num2str(begin_of_signal),' mel']);
    disp(['end: ',num2str(end_of_signal),' mel']);
  elseif filetype==8
    disp(['begin: ',num2str(begin_of_signal),' bark']);
    disp(['end: ',num2str(end_of_signal),' bark']);
  elseif filetype==9
    disp(['begin: ',num2str(begin_of_signal),' ERB']);
    disp(['end: ',num2str(end_of_signal),' ERB']);
  elseif filetype==10
    disp(['begin: ',num2str(begin_of_signal),' log10( f / Hz )']);
    disp(['end: ',num2str(end_of_signal),' log10( f / Hz)']);
  end;
  if n_lines==1
    error('selected range too small');
  else
     disp(['number_of_lines: ', num2str(n_lines)]);
  end;
  
% cut matrices and xaxis-vector
  if isempty(xaxis)==0
    xaxis = xaxis(:,first_line:last_line);
  end;
  if isempty(first_channel)==0
    first_channel = first_channel(:,first_line:last_line);
  end;
  if isempty(second_channel)==0
    second_channel = second_channel(:,first_line:last_line);
  end;
  if isempty(first_phase)==0
    first_phase = first_phase(:,first_line:last_line);
  end;
  if isempty(second_phase)==0
    second_phase = second_phase(:,first_line:last_line);
  end;
end;
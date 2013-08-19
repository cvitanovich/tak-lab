function answer=askyesno(prompt);
% ASKYESNO(prompt) gets a 'yes' (returns 1) or 'no'
%         (returns 0) answer from the user. Only the
%         first character of the input is recognized.

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 26.08.1994 A. Raab

while 1
  answer=lower(input(prompt,'s'));
  if ~isempty(answer)
    answer=find(answer(1)==['ny'])-1;
  end;
  if ~isempty(answer) break, end;
end;

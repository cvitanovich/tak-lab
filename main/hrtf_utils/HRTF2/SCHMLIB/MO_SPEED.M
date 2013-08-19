function fps = mo_speed()

% function fps = mo_speed()
%
% MO_SPEED changes the speed of motions. The user input must be an positive
%          number, and defines how many frames per second will be displayed.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

fps = input('How many figures per second ? ');
if isempty(fps)
  fps = 1;
end;

if fps<0
  fps = -fps;
end;
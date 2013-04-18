function [h] = findax(fig)

% function [h] = findax(fig)
%
% FINDAX returns a vector 'h' with the handles of all axes, that are
%        situated in the figure with the handle 'fig'.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

h = [];
% find axes objects
c = get(fig,'Children');
j=1;         %axis counter
for i = 1:length(c),
    if strcmp(get(c(i),'Type'),'axes')
      h(j) = c(i);
      j = j + 1;
    end;
end;
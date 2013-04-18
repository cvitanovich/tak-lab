function scaling = frqscale()

% function scaling = frqscale()
%
% SCALING returns a string 'scaling' defining the kind of scaling for the
%         frequency axis. The scaling can be linear ('lin') or logarithm
%         ('log'), what the user has to decide. Explicit scaling is used for
%         contour and mesh plots, where separate plot functions are not provided.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

scaling = [];
while isempty(scaling)
  scaling = input('''lin'' or ''log'' scaling of the frequency axis ? ','s');
  [row,letters] = size(scaling);
  if letters~=3
    scaling = [];
  elseif ( (scaling~='lin')&(scaling~='log') )
   scaling = [];
  end;
end;
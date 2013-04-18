function s=dinput(prompt,default,mode);
% DINPUT Keyboard input with default value.
%        s=DINPUT(prompt,default,mode)
%        Invokes the MATLAB function INPUT with parameters
%        prompt and mode. If the returned value from INPUT 
%        is empty DINPUT returns the value default.

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 17.08.1994 A. Raab

if nargin<3 
  s=input(prompt);
else
  s=input(prompt,mode);
end;
if isempty(s) s=default; end;ÿ

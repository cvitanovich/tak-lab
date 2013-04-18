function sel=textmenu(header,arg1,arg2,arg3,arg4,arg5,arg6,...
                             arg7,arg8,arg9);
% TEXTMEU display a text menu and returns the user selection.
%         sel=TEXTMENU(header,arg1,arg2,arg3,arg4,arg5,arg6,...
%                             arg7,arg8,arg9);

disp(header);
disp(' ');
for i=1:nargin-1
  eval(sprintf('s=arg%i;',i));
  disp(sprintf('  %i. %s',i,s));
end;
sel=dinput('Enter command: ',0);

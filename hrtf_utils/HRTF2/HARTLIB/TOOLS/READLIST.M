function list=readlist(cmdfile)
% READLIST reads a ascii file containing a file list
%          list=READLIST(cmdfile)

% Achim Raab
% Lehrstuhl f. allg. Elektrotechnik und Akustik
% Ruhr-Universtitaet Bochum
% 1994
% Modified
% 12.5.1996 Klaus  Hartung

disp('Reading List. May take a while ...')
fid=fopen(cmdfile,'r');

while 1
  line=fgetl(fid);
%  disp(line);
  if line(1)==-1 break; end;
  if ~isempty(line)
    line=fliplr(deblank(fliplr(line)));
    list=str2mat(list,line);
  end;
end;
fclose(fid);
list=list(2:size(list,1),:);

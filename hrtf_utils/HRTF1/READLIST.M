function list=readlist(cmdfile)
% READLIST reads a ascii file containing a file list
%          list=READLIST(cmdfile)
fid=fopen(cmdfile,'r');

list =[];

while 1
  line=fgetl(fid);
  disp(line);
  if line(1)==-1 break; end;
  if ~isempty(line)
    line=fliplr(deblank(fliplr(line)));
    list=str2mat(list,line);
  end;
end;
fclose(fid);
list=list(2:size(list,1),:);

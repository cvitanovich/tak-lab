function ch=seldirs(dirs,ch);
% SELDIRS displays the available directions (44!)
%         and enables the user to select special
%         directions.
%         ch=SELDIRS(dirs,presel);

ndirs=size(dirs,1);

if nargin<2
  ch=ones(size(dirs,1),1);
end;

chr=['-x'];
xm=[1:10;11:20;21:30;31:40];
flag=0;

if ndirs~=44 error('number of directions must be 44'); end;

while flag==0
  clc;
  for i=1:size(xm,1)
    s=sprintf('az %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i',dirs(xm(i,:),1));
    disp(s);
    s=sprintf('el %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i',dirs(xm(i,:),2));
    disp(s);
    s=sprintf('   %4c %4c %4c %4c %4c %4c %4c %4c %4c %4c %4c',chr(ch(xm(i,:))+1));
    disp(s);
  end;
  s=sprintf('az %4i %4i %4i %4i',dirs(41:44,1));
  disp(s);
  s=sprintf('el %4i %4i %4i %4i',dirs(41:44,2));
  disp(s);
  s=sprintf('   %4c %4c %4c %4c',chr(ch(41:44)+1));
  disp(s);

  c=dinput('Enter index of directions to toggle, 0 to exit [1,3:21]: ',[1,3:21]);
  if c==0 
    flag=1;
  else
    if max(c)<=ndirs & min(c)>0  ch(c)=~ch(c); end;
  end;
end;


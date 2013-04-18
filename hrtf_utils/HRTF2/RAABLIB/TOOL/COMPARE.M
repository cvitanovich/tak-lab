function compare(catname1,catname2)

N=2048;

s=['load ' catname1];
eval(s);
B1=B;
if exist('A')~=1
  A=ones(size(B,1),1);
end;
A1=A;

s=['load ' catname2];
eval(s);
B2=B;
if exist('A')~=1
  A=ones(size(B,1),1);
end;
A2=A;

for k=1:size(B1,1)
  j=round(k/2);
  H1=freqz(B1(k,:),A1(k,:),N);
  H2=freqz(B2(k,:),A2(k,:),N);
  gtitle=sprintf('az=%i, el=%i, channel=%i',dirs(j,1),dirs(j,2),iseven(k));
  magplot([H1(:) H2(:)],44100,gtitle,'rg',0,15000);
  grid on;
  grcntrl;
  legend('r',catname1,'g',catname2,-1);
  if ~askyesno('Continue (y/n)') break; end;
  close;
end;
function cof2sr2(coffile,sr2nameprefix)
% COF2SR2 converts a coffile with filter coefficients
%         to a series of SR2 files. The name of each SR2 file
%         is starting with the given string sr2nameprefix
%         followed by a number for the direction.
%         Each processed direction is printed on the screen.
%         COF2SR2(coffile,sr2nameprefix)

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 15.10.1994 A. Raab

eval(['load ' coffile]);
K=size(B,1);
L=size(B,2);
if exist('A')==1              % its a IIR filter
  delta=zeros(1,1024);        % compute impulse response
  delta(1)=1;
  H=zeros(K,1024);
  for i=1:K
    H(i,:)=filter(B(i,:),A(i,:),delta);
  end;
else
  H=[B zeros(K,1024-L)];
end;

j=1;
disp('No.    az     el   file');
for i=1:2:K
  fname=sprintf('%s_%i.sr2',sr2nameprefix,j);
  s=sprintf('%2i    %3i    %3i   %s',j,dirs(j,1),dirs(j,2),fname);
  disp(s);
  writesr2(fname,H(i,:),H(i+1,:));
  j=j+1;
end;

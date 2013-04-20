function gridhandle=usergrid(X,Y,C);
% USERGRID      Draw grid lines 
%               USERGRID(X,Y,C) draws dashed lines at specified 
%               positions with color C and returns a handle to all
%               lines. Y and C are optional.

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 03.07.1994 A. Raab

if nargin<3 C='w'; end;
if nargin<2 Y=[]; end;

gridhandle=zeros(1,length(X)+length(Y));
ax=axis;
for i=1:length(X)
  gridhandle(i)=line([X(i) X(i)],[ax(3) ax(4)]);
  set(gridhandle(i),'LineStyle','--','Color',C);
end;
for j=1:length(Y)
  gridhandle(j+i)=line([ax(1) ax(2)],[Y(j) Y(j)]);
  set(gridhandle(j+i),'LineStyle','--','Color',C);
end;

function semlogx2(x,y,s)

if (nargin<3) s='y'; end

maxgrid=150;
nticks=9;


stdoct=[100 200 500 1000 2000 5000 10000 20000 50000];

%minstdoct=min(stdoct);
%maxstdoct=max(stdoct);

maxx=max(x);
minx=min(x);
if (minx == 0) 
 minx=x(2); 
 n_grid=maxgrid-1;
 else
 n_grid=maxgrid;
end;

%if (minx == 0) minx=x(2); end;
lx=length(x);

k=find(stdoct<=maxx);
stdoct=stdoct(k);
k=find(stdoct>=minx);
stdoct=stdoct(k);
nticks=length(stdoct);

maxoct=log(maxx/minx)/log(2);
xlog=(2.^linspace(0,maxoct,n_grid))*minx;

if (x(1)==0) xlog=[0 xlog]; end;

%xlog=(2.^linspace(0,maxoct,maxgrid-1))*minx;

%xlog=[minx(1) xlog];



xticklabels=zeros(nticks,5);

for nt=1:nticks
  %xtick(nt)=maxgrid/nticks*nt; 
  %label=int2str(fix(maxx/2^(maxoct/nticks*(nt-1))));
  %sizelabel=length(label);
  %xticklabels(nticks-nt+1,1:sizelabel)=label;
 
  diff=xlog-stdoct(nt)*ones(size(xlog));
  [v,p]=min(abs(diff));
  xtick(nt)=p;
  label=int2str(fix(stdoct(nt)));
  sizelabel=length(label);
  xticklabels(nt,1:sizelabel)=label;  
  
end;

yi=interp1(x,y,xlog);

px=1:maxgrid;
plot(px,yi,s);
set(gca,'XTick',xtick);
set(gca,'XTickLabels',xticklabels);


drawnow;

set(gcf,'UserData',[ maxgrid xlog minx maxx yi]);


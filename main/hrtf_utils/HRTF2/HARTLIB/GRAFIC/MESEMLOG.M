function mesemlog(x,y,z)

%maxgrid=400;
maxgrid=150
nticks=9;


%stdoct=[120 240 480 960 2000 4000 8000 16000];
stdoct=[100 200 500 1000 2000 5000 10000 20000 50000];

maxx=max(x);
minx=min(x);
if (minx == 0) 
 minx=x(2); 
 n_grid=maxgrid-1;
 else
 n_grid=maxgrid;
end;

lx=length(x);
dx=x(2)-x(1);

k=find(stdoct<=maxx);
stdoct=stdoct(k);
k=find(stdoct>=minx);
stdoct=stdoct(k);
nticks=length(stdoct);

maxoct=log(maxx/minx)/log(2);
xlog=round((2.^linspace(0,maxoct,n_grid))*minx);
x=round(x);
if (x(1)==0) xlog=[0 xlog]; end

%colormap('jet')

set(gcf,'UserData',[ maxgrid xlog minx maxx]);


set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperType','a4letter');
p=get(gcf,'PaperPosition');
%p = [2.5    2.5   26   16];
p(1)=2.5;
p(2)=2.5;
p(3)=16;
p(4)=14;
set(gcf,'PaperPosition',p);

xticklabels=zeros(nticks,5);

for nt=1:nticks
 % xtick(nt)=maxgrid/nticks*nt; 
 % label=int2str(fix(maxx/2^(maxoct/nticks*(nt-1))));
 % sizelabel=length(label);
 % xticklabels(nticks-nt+1,1:sizelabel)=label;

  diff=xlog-stdoct(nt)*ones(size(xlog));
  [v,p]=min(abs(diff));
  xtick(nt)=p;
  label=int2str(fix(stdoct(nt)));
  sizelabel=length(label);
  xticklabels(nt,1:sizelabel)=label;  
  
end;

[mz,nz]=size(z);

%max(x)
%min(x)

%max(xlog)
%min(xlog)
for cz=1:nz
   zi(:,cz)=interp1(x,z(:,cz),xlog);
end


%zi=interp1(x,z,xlog);

px=1:maxgrid;
%waterfall(px,y,zi');
%mesh(px,y,zi',zeros(size(zi')));
%size(px)
%size(y)
%size(xlog)
%size(zi')

%mesh(px,y,zi');

py=min(y):5:max(y);
[miz,niz]=size(zi);

for lz=1:miz
  zi2(lz,:)=(interp1(y,zi(lz,:)',py)).';
end

%size(px)
%size(py)
%size(zi2)
mesh(px,py,zi2');

%waterfall(px,y,zi');

%meshz(px,y,zi');
set(gca,'FontSize',16);

set(gca,'XTick',xtick);
set(gca,'XTickLabels',xticklabels);
drawnow;
%set(gcf,'UserData',[ maxgrid xlog minx maxx zi px]);

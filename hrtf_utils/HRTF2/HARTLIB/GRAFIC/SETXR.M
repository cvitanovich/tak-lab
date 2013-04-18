
function setxr(xmin,xmax)

ui_handle=get(gcf,'UserData');

maxgrid=ui_handle(1);
xlog=ui_handle(2:2+maxgrid-1);
minx=ui_handle(2+maxgrid);
maxx=ui_handle(2+maxgrid+1);
%zi=ui_handle(2+maxgrid+2:2+maxgrid+2+maxgrid-1);
%px=ui_handle(4+2*maxgrid:4+2*maxgrid+maxgrid-1);

diff=xlog-xmin;
[k,p]=min(abs(diff));
x(1)=p;
diff=xlog-xmax;
[k,p]=min(abs(diff));
x(2)=p;

set(gca,'XLim',x);


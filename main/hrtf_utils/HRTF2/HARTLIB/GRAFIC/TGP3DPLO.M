function tg3dplo(H,Fs,minf,maxf,y);
docorrect=0;

[a,b]=size(H);

if (nargin < 3) minf=0; end;
if (nargin< 4) maxf=Fs; end
if (nargin < 5) y=1:b; end;

f=linspace(0,Fs/2,length(H(:,1)))

[v,k]=min(abs(f-minf*ones(size(f))));
minx=k;

[v,k]=min(abs(f-maxf*ones(size(f))));
maxx=k;


for k=1:b
 TG(:,k)=tgroup(H(:,k));
end

mesemlog(f(minx:maxx),y,TG(minx:maxx,:))
xlabel('frequency/Hz');
zlabel('phase delay/s');
%set(gca,'ZLim',[min(min(TG(minx:max,:))),max(max(TG(minx:maxx,:)))]);
set(gca,'YLim',[min(y), max(y)]);
view(45,30);

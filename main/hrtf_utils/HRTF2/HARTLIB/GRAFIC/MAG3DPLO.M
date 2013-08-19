function mag3dplo(H,Fs,minf,maxf,y);

[a,b]=size(H);

if (nargin < 3) minf=0; end;
if (nargin< 4) maxf=Fs; end
if (nargin < 5) y=1:b; end;

f=linspace(0,Fs/2,length(H(:,1)))

[v,k]=min(abs(f-minf*ones(size(f))));
minx=k

[v,k]=min(abs(f-maxf*ones(size(f))));
maxx=k
f(maxx)
Hlog=20*log10(abs(H(minx:maxx,:)));
mesemlog(f(minx:maxx),y,Hlog)
xlabel('frequency/Hz');
zlabel('magnitude/dB');
set(gca,'ZLim',[floor(min(min(Hlog))/10)*10,ceil(max(max(Hlog))/10)*10]);
%set(gca,'ZLim',-60,20);
set(gca,'YLim',[min(y), max(y)]);
%h=colorbar('vert');
%t=get(h,'title');
%set(t,'string','dB');
grid on
view(45,45);

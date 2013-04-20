function meantf(filelist,outname,t)
% MEANTF computes the mean magnitude transfer function for the 
%        specified file list and writes it to the file outname.
%        MEANTF(filelist,outname)
%
epsilon=1E-8;

if (nargin < 3) t = 2.5 ;end;

samplerate=24000;
fftlen=1024;
chname='lr';
epsilon=1E-8;
n=t*samplerate/1000.0;

Mmean=zeros(fftlen,2);          % mean magnitude both channels
list=readlist(filelist);
nfiles=size(list,1);

disp('computing mean magnitude function...');
for i=1:nfiles
  filename=list(i,:);
  [Hl,Hr]=readsc2(filename);
   
  Hl=Hl(:);
  Hr=Hr(:);


  [hl,hr] = luxusfft(Hl,Hr,n);

  Hl = rfft(hl,2048);
  Hr = rfft(hr,2048);

  Mmean(:,1)=Mmean(:,1)+abs(Hl);
  Mmean(:,2)=Mmean(:,2)+abs(Hr);
end;

Mmean=Mmean/nfiles;
info=['Mean magnitude transfer function for ' filelist];
eval(['save ' outname ' Mmean info;']);

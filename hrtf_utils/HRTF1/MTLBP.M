function mtlbp(outfile,f1,f2,sampling_rate,n_lines)

% function mtlbp(outfile,f1, f2, sampling_rate,n_lines)

if(nargin <3) error('Args are missing'); return; end;
if(nargin <4) sampling_rate=44100.0; end
if(nargin <5) n_lines=1024; end;

n_lines=n_lines*2;

filetype=3;
info_blocks=2;
n_channels=1;
first_line=0;
resolution=sampling_rate/(n_lines);
last_line=resolution*(n_lines-1);
n_directions=1;
comment1='Bandpass Filter';
comment2=[];
mtlwh(outfile,filetype,info_blocks,n_channels,n_lines/2, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);

f_0=f1/3.
f_3=1.1*f2
res=sampling_rate/n_lines
l_0=fix(f_0/res)+1;
l_3=fix(f_3/res)+1;

b=fir1(100,[f1/(sampling_rate/2) f2/(sampling_rate/2)],blackman(100+1));
irplot(b',sampling_rate);
B=rfft(b,n_lines);
l_0
l_3
%H(1:n_lines/2)=ones(1,n_lines/2);
B(1,1:l_0)=ones(1,l_0)*1E-18;
B(1,l_3:(n_lines/2))=ones(1,(n_lines/2)-l_3+1)*1E-18;
n_lines
size(B)
%size(H)

B=B.*H;
magplot(B',sampling_rate,'','',100,12000);
phasplot(B',sampling_rate,'','',100,12000);

mtlwch(outfile,B,1);

return;

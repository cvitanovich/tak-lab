%function []=cleanhrir(filename,t);

% cleanhrir cleans the measured impulse response   
%          CLEANHRIR(filename,t);
%          t estimated length of impulse response
%          PRELIMINARY VERSION !!

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
%     RUB
% (p) A. Raab, K. Hartung, 09.12.1994, 8.00 pm
%  V 1.2 K. Hartung 13.12.94 2.30 pm

%if (nargin < 2) t = 2.5 ;end;

[filetype,info_blocks,...
          n_channels,n_lines, sampling_rate,...
          first_line,last_line,n_directions,...
          comment1, comment2] = readh(filename)

sampling_rate=24000;
epsilon=1E-8;
n=t*samplerate/1000.0;

[Hl,Hr,azi,ele]=readsc2(filename);
Hl=Hl(:);
Hr=Hr(:);

load 802B_DTF
[hl,hr] = luxusfft(Hl,Hr,n,Mmean);

%Hl = rfft(hl,2048);
%Hr = rfft(hr,2048);

%Hinteraural=Hl(:)./Hr(:);
%Hinteraural(1:171)=abs(real(Hinteraural(1:171)));

strrep(filename,'.SC2','.sr2')
filetype=2;
n_lines=length(hl);
first_line=0;
last_line=(n_lines-1)/sampling_rate;

writeh(newfilename,filetype,info_blocks,...
          n_channels,n_lines, sampling_rate,...
          first_line,last_line,n_directions,...
          comment1, comment2)
wrt_d(newfilename,[ele;azi]);

[fid,msg]=fopen(newfilename,'a','ieee-le');
if fid<0 error(msg); end;
fseek(fid,256*info_blocks,'bof');
fwrite(fid,hl,'float32');
fwrite(fid,hr,'float32');
fclose(fid);

filename='m1.sc2'

[filetype, info_blocks,...
n_channels, n_lines, sampling_rate,...
first_line, last_line, n_directions,...
comment1, comment2] = mtlrh(filename);

filesize=sizeoffile(filename);

n_lines_v = (filesize-256*info_blocks)/(n_channels*4*2);
N_LINES=n_lines_v;


for n=1:11

filename=[ 'm' int2str(n) '.sc2'];

HRTF_RAW(:,n*2-1)=mtlrch(filename,1);
HRTF_RAW(:,n*2)=mtlrch(filename,2);

end


HRIR_RAW = rifft(HRTF_RAW);
HRIR_RAW = HRIR_RAW/max(max(HRIR_RAW))*32767;

HRIR     = optwin(HRIR_RAW,128,512);

HRTF     = rfft(HRIR);

IATF = HRTF(:,1:2:2*n)./HRTF(:,2:2:2*n);



MAP = [-10:10:90; zeros(1,11)]



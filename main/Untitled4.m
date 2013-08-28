ATTENS=1.*ones(1,100);
SCALES=2.*ones(1,100);
VOLTAGES=3.*ones(1,100);
fid = fopen('test3.txt','w');
DATA=[ATTENS; SCALES; VOLTAGES];
formatSpec='%3.1f\t%5.0f\t%5.3f\n';
nbytes=fprintf(fid,formatSpec,DATA);
fclose(fid);
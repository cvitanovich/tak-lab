% load knowles header file
uiopen('~/data/calib/*.mat');

% load raw data (Freefield)
FilterSpec='/home/andrew/data/calib/*.txt';
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec);
F=dlmread([PathName FileName],'\t'); % read tab delimited data to matrix

% load raw data (Intraural)
FilterSpec='/home/andrew/data/calib/*.txt';
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec);
I=dlmread([PathName FileName],'\t'); % read tab delimited data to matrix


attens=unique(F(:,1));
% remove data with really small rms values
[rowsF,cols]=find(F(:,3)>=KNOWLES.coefs(3));
[rowsI,cols]=find(I(:,3)>=KNOWLES.coefs(3));
rows=intersect(rowsF,rowsI);
F=F(rows,:);
dBs = (1/KNOWLES.coefs(2))*log( (F(:,3) - KNOWLES.coefs(3)) ./ KNOWLES.coefs(1) );
F=[F dBs];
I=I(rows,:);
dBs = (1/KNOWLES.coefs(2))*log( (I(:,3) - KNOWLES.coefs(3)) ./ KNOWLES.coefs(1) );
I=[I dBs];

% remove all measurements with dB > 55 or < 30
upper=55; lower=30;

[rowHi,cols]=find(F(:,4)<upper);
[rowLo,cols]=find(F(:,4)>lower);
rowF=intersect(rowHi,rowLo);
[rowHi,cols]=find(I(:,4)<upper);
[rowLo,cols]=find(I(:,4)>lower);
rowI=intersect(rowHi,rowLo);
rows=intersect(rowF,rowI);
GAIN=I(rows,4)-F(rows,4);
mu=mean(GAIN);
stdev=std(GAIN);
str=['Mean = ' num2str(mu) ' STDEV= ' num2str(stdev)];
figure; boxplot(GAIN); hold on;
title(str);
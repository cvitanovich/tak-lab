clear all
Fs = 30000;

[r1000 phi1000] = lonoisenoise(1000,Fs);
save carrier1k.mat

clear all
Fs = 30000;
s = RandStream('mt19937ar','Seed',1);
[r1000 phi1000] = lonoisenoiseNorm(1000,Fs,s);
save carrier1k_norm1.mat

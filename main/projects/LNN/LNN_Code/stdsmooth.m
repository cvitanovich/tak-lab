clear all
Fs = 30000;
s = RandStream('mt19937ar','Seed',1 + (2^32 - 1).*rand(1,1));
[r1000 phi1000] = LNNstd(1000,Fs,s);
save carrier1k_smooth1.mat

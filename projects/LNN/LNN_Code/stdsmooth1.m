clear all
Fs = 30000;
a=clock;
s = RandStream('mt19937ar','Seed',round(a(end)*10^5));
[s500 phi500] = LNNstd(500,Fs,s);
save carrier1k_smooth1.mat

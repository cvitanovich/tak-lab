function [] = plot_stims(L,R);

%function [] = plot_stims(L,R);

global FN

Fs = 30000;							% sampling frequency
cF = round(1000*exp(([12:40]/12)*log(2)))';
weights = ones(size(cF));

tempL = ERBFilterBankA(L, cF, Fs);		% has dimensions n_cF x length(noi)
tempR = ERBFilterBankA(R, cF, Fs);

% calc and plot ILD
ILDnew = calclevel_time(tempL, tempR, cF);
figure(6); hold on;
plot(cF,ILDnew,'r');

% calc and plot ITD
[ITDnew, time] = calcitd(tempL,tempR, cF, Fs, weights);
figure(7); hold on;
plot(time,ITDnew,'r');

% add to saved file
if exist(['''' FN.temp_stim_path  'XXX.mat' ''''])
    eval(['load ' FN.temp_stim_path 'XXX'])
    [m,n] = size(ILD);
    ILD(m+1,:) = ILDnew;
    [m,n] = size(ITD);
    ITD(m+1,:) = ITDnew;
else
    ITD = ITDnew;
    ILD = ILDnew;
end
eval(['save ' FN.temp_stim_path 'XXX ILD ITD'])

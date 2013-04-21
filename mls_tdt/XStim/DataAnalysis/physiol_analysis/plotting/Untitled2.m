position = [300 0 950 950];
xvals = 2000:450:11000;
yvals = -30:4:30;
xtext = 'Freq (Hz)';
ytext = 'ILD (dB)';
xoffset = 5;
yoffset = 5;
numreps = 10;
tot_dur = 200;
stim_dur = 100;
startcount = 0;
endcount = 100;
row_space = 5;
col_space = .2*tot_dur;

spikedata = round(150 * rand(1,100));


[hrasterfig,hrastaxes] = makerasterfig(position,...
   xvals,...
   yvals,...
   xtext,...
   ytext,...
   xoffset,...
   yoffset,...
   numreps,...
   tot_dur,...
   stim_dur,...
   startcount,...
   endcount,...
   row_space,...
   col_space);

plotraster(hrasterfig,spikedata,2,3,5,tot_dur,stim_dur,numreps,xoffset,yoffset,row_space,col_space);
hold
plotraster(hrasterfig,spikedata,2,3,4,tot_dur,stim_dur,numreps,xoffset,yoffset,row_space,col_space);

%% INITIATE TDT PARAMETERS HERE
TDT.nPlayChannels=2;
TDT.nBuffers=2;
TDT.bufpts=[npts npts];
TDT.nRecChannels=2;
TDT.dec_factor=5;
TDT.din = 1;
TDT.Fs = 30000;
TDT.npts_total_play=???;

%% INITIALIZE TDT
out=TDT_init;
if out<0
	return;
end

%% INITIALIZE BUFFERS
TDT_buffers;

%% INITIALIZE PD1
PD1_init;

%% SET ATTENS
TDT_attens(attens);




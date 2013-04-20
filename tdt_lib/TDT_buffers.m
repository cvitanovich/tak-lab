function TDT=TDT_buffers(TDT)
% generic setup for stimulus play specifications for TDT Sys II
% permits nearly any combination of of play/record channel setups
% and an unlimited number of buffers in the play sequence
% with variable combinations of buffer lengths for each channel

% input (a structure):
% TDT.nPlayChannels -- # of channels desired (2 for left/right headphones, etc.)

% TDT.playpts -- a cell array of bufpts for each buffer in each channel,
% formatted as such: { CHANNEL 1, CHANNEL 2, ...}
% each cell contains the pts in the buffer.
% for example: CHANNEL 1 = [ch1_bufpts1 ch1_bufpts2 ...]

% TDT.nRecChannels -- # recording channels desired

% TDT.recpts -- a cell array of bufpts for each recording buffer in each
% channel, formatted just like TDT.bufpts

% TDT.dec_factor -- decimation factor


% outputs:
% TDT.play_spec -- buffer id for play specification
% TDT.ch_seq -- buffer id for each channel play sequence

% TDT.stim_buffers - buffer id values for each buffer
% formatted as cell array just like TDT.playpts

% TDT.rec_spec -- recording spec buffer
% TDT.rec_ch_seq -- recording channel sequence

% TDT.rec_buffers -- id values for each recording buffer
% formatted just like TDT.recpts

% TDT.dec_buffers -- id values for decimated recording buffers
% formatted just like TDT.recpts

% TDT.n_total_buffers -- total no. of buffers in use

% PLAY SEQUENCE SETUP
TDT.play_spec=1; S232('allot16',TDT.play_spec,(TDT.nPlayChannels+1));

% ALLOT PLAY CHANNEL SEQUENCE BUFFERS
TDT.ch_seq=2:(TDT.nPlayChannels+1);
for ch=1:length(TDT.ch_seq)
	S232('allot16',TDT.ch_seq(ch),(TDT.nPlayChannels+1));
end

% ALLOT PLAY BUFFERS
TDT.stim_buffers=cell(1,TDT.nPlayChannels);
nbuffers=length(TDT.playpts{1});
start=(TDT.ch_seq(end)+1); stop=(TDT.ch_seq(end)+nbuffers);
for ch=1:TDT.nPlayChannels
    
	TDT.stim_buffers{ch}=start:stop;
	for buf=1:nbuffers
		S232('allot16',TDT.stim_buffers{ch}(buf),TDT.playpts{ch}(buf));
    end
    
    if ch<TDT.nPlayChannels
        nbuffers=length(TDT.playpts{ch+1});
        start=stop+1;
        stop=stop+nbuffers;
    end
end

% PLAY SPECIFICATION LIST
S232('dpush',(TDT.nPlayChannels+1));
S232('value',0);
for j=0:(TDT.nPlayChannels-1)
	S232('make',j,TDT.ch_seq(j+1));
end
S232('make',TDT.nPlayChannels,0);
S232('qpop16',TDT.play_spec);

% PLAY SEQUENCES
for ch=1:TDT.nPlayChannels
    nbuffers=length(TDT.playpts{ch});
	S232('dpush',(nbuffers*2+1));
	S232('value',0);
	cnt=0;
	for buf=1:nbuffers
		S232('make',cnt,TDT.stim_buffers{ch}(buf));
		cnt=cnt+1;
		S232('make',cnt,1);
		cnt=cnt+1;
	end
	S232('make',cnt,0);
	S232('qpop16',TDT.ch_seq(ch));
end

% RECORDING SETUP
if(TDT.nRecChannels)
	% recording specification
	N=TDT.stim_buffers{TDT.nPlayChannels}(end);
	TDT.rec_spec=N+1; S232('allot16',TDT.rec_spec,(TDT.nRecChannels+1));

	% recording channels
    N=N+2;
	TDT.rec_ch_seq=N:(N + TDT.nRecChannels - 1);
	for ch=1:length(TDT.rec_ch_seq)
		S232('allot16',TDT.rec_ch_seq(ch),(TDT.nRecChannels+N+1));
	end

	% recording buffer sequence
	N=TDT.rec_ch_seq(end);
	TDT.rec_buffers=cell(1,TDT.nRecChannels);
    nbuffers=length(TDT.recpts{1});
	start=(N+1); stop=(N+nbuffers);
	for ch=1:TDT.nRecChannels
        TDT.rec_buffers{ch}=start:stop;
		for buf=1:nbuffers
			S232('allot16',TDT.rec_buffers{ch}(buf),TDT.recpts{ch}(buf));
        end
        if(ch<TDT.nRecChannels)
            nbuffers=length(TDT.recpts{ch+1});
            start=stop+1;
            stop=stop+nbuffers;
        end
    end
    
    
    
    
    
    if(TDT.dec_factor) % decimating?
        % decimated record buffers
        N=TDT.rec_buffers{TDT.nRecChannels}(end);
        TDT.dec_buffers=cell(size(TDT.rec_buffers));
        nbuffers=length(TDT.recpts{1});
        start=(N+1); stop=(N+nbuffers);
        for ch=1:TDT.nRecChannels
            TDT.dec_buffers{ch}=start:stop;
            for buf=1:nbuffers
                dec_pts = ceil(TDT.recpts{ch}(buf) / 2^TDT.dec_factor);
                S232('allot16',TDT.dec_buffers{ch}(buf),dec_pts);
            end
            if(ch<TDT.nRecChannels)
                nbuffers=length(TDT.recpts{ch+1});
                start=stop+1;
                stop=stop+nbuffers;
            end
        end
        N=TDT.dec_buffers{TDT.nRecChannels}(end);
    end
    
	% record spec list
	S232('dpush',(TDT.nRecChannels+1));
	S232('value',0);
	for j=0:(TDT.nRecChannels-1)
		S232('make',j,TDT.rec_ch_seq(j+1));
	end
	S232('make',TDT.nRecChannels,0);
	S232('qpop16',TDT.rec_spec);

	% record sequences for each recording channel
	for ch=1:TDT.nRecChannels
        nbuffers=length(TDT.recpts{ch});
		S232('dpush',(nbuffers*2+1));
		S232('value',0);
		cnt=0;
		for buf=1:nbuffers
			S232('make',cnt,TDT.rec_buffers{ch}(buf));
			cnt=cnt+1;
			S232('make',cnt,1);
			cnt=cnt+1;
		end
		S232('make',cnt,0);
		S232('qpop16',TDT.rec_ch_seq(ch));
	end
end

TDT.n_total_buffers=N;


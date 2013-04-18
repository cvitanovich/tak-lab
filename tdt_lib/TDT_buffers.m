function TDT_buffers
% generic setup for stimulus play specifications for TDT Sys II
% permits any number of play/record channel combinations
% and an unlimited number of buffers in the play sequence
global TDT

% input (a structure):
% TDT.nPlayChannels -- # of channels desired (2 for left/right headphones, etc.)
% TDT.nBuffers -- # of stimulus buffers in the play sequence
% TDT.bufpts -- no. pts in each buffer (1 x TDT.nBuffers vector)
% TDT.nRecChannels -- # recording channels desired
% TDT.dec_factor -- decimation factor

%
% outputs:
% TDT.play_spec -- buffer id for play specification
% TDT.ch_seq -- buffer id for each channel play sequence
% TDT.stim_buffers - buffer id values for each buffer
% TDT.rec_spec -- recording spec buffer
% TDT.rec_ch_seq -- recording channel sequence
% TDT.rec_buffers -- id values for each recording buffer
% TDT.dec_buffers -- id values for decimated recording buffers
% TDT.n_total_buffers -- total no. of buffers in use

% PLAY SEQUENCE SETUP
TDT.play_spec=1; S232('allot16',TDT.play_spec,(TDT.nPlayChannels+1));

% ALLOT PLAY CHANNEL SEQUENCE BUFFERS
TDT.ch_seq=2:(TDT.nPlayChannels+1);
for ch=1:length(TDT.ch_seq)
	S232('allot16',TDT.ch_seq(ch),(TDT.nPlayChannels+1));
end

% ALLOT PLAY BUFFERS
TDT.stim_buffers=NaN*ones(TDT.nPlayChannels,TDT.nBuffers);
start=(TDT.ch_seq(end)+1); stop=(TDT.ch_seq(end)+TDT.nBuffers);
for ch=1:TDT.nPlayChannels
	TDT.stim_buffers(ch,:)=start:stop;
	start=stop+1;
    for j=1:TDT.nBuffers
        stop=stop+j;
    end
	for buf=1:TDT.nBuffers
		S232('allot16',TDT.stim_buffers(ch,buf),TDT.bufpts(buf));
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
	S232('dpush',(TDT.nBuffers*2+1));
	S232('value',0);
	cnt=0;
	for buf=1:TDT.nBuffers
		S232('make',cnt,TDT.stim_buffers(ch,buf));
		cnt=cnt+1;
		S232('make',cnt,1);
		cnt=cnt+1;
	end
	S232('make',cnt,0);
	S232('qpop16',TDT.ch_seq(ch));
end

% RECORDING SETUP
if TDT.nRecChannels
	% recording specification
	N=TDT.stim_buffers(end);
	TDT.rec_spec=N+1; S232('allot16',TDT.rec_spec,(TDT.nRecChannels+1));

	% recording channels
    N=TDT.rec_spec(end);
	TDT.rec_ch_seq=(N+1):(TDT.nRecChannels+N+1);
	for ch=1:length(TDT.rec_ch_seq)
		S232('allot16',TDT.rec_ch_seq(ch),(TDT.nRecChannels+N+1));
	end

	% recording buffer sequence
	N=TDT.rec_ch_seq(end);
	TDT.rec_buffers=NaN*ones(TDT.nRecChannels,TDT.nBuffers);
	start=(N+1); stop=(N+TDT.nBuffers);
	for ch=1:TDT.nRecChannels
		TDT.rec_buffers(ch,:)=start:stop;
		start=stop+1;
		stop=stop+TDT.nBuffers;
		for buf=1:TDT.nBuffers
			S232('allot16',TDT.rec_buffers(ch,buf),TDT.bufpts(buf));
		end
	end
    
    N=TDT.rec_buffers(end);
    if(TDT.dec_factor) % decimating?
        % decimated record buffers
        
        start=(N+1); stop=(N+TDT.nBuffers);
        for ch=1:TDT.nRecChannels
            TDT.dec_buffers(ch,:)=start:stop;
            start=stop+1;
            stop=start+TDT.nBuffers;
            for buf=1:TDT.nBuffers
                dec_pts = ceil(TDT.bufpts(buf) / 2^TDT.dec_factor);
                S232('allot16',TDT.dec_buffers(ch,buf),dec_pts);
            end
        end
        N=TDT.dec_buffers(end);
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
		S232('dpush',(TDT.nBuffers*2+1));
		S232('value',0);
		cnt=0;
		for buf=1:TDT.nBuffers
			S232('make',cnt,TDT.rec_buffers(buf));
			cnt=cnt+1;
			S232('make',cnt,1);
			cnt=cnt+1;
		end
		S232('make',cnt,0);
		S232('qpop16',TDT.rec_ch_seq(ch));
	end
end
TDT.n_total_buffers=N;


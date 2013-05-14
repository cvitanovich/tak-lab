%% INIT TDT
status=0;
while ~status
	S2init=S232('S2init', 0, 'INIT_SECONDARY', 20000);
	APlock=S232('APlock',200, 0);
	XBlock=S232('XBlock',200, 0);
	if S2init && APlock && XBlock
		status=1;
	else
		s232('APunlock', 0);
        s232('S2close');
        err = S232('getS2err');
        switch err
            case 0
            case 1
                h=warndlg('APOS error in initiation','warning'); 
                uiwait(h);
            case 2
                h=warndlg('XBUS error in initiation','warning'); u
                iwait(h);
        end
		err='Cannot initiate TDT. Try resetting, then click okay.';
		choice=questdlg(err,'S232 Error','OKAY','ABORT','OKAY');
		if strcmp(choice,'ABORT')
			out=-1; return;
		end
	end
end
S232('trash');
S232('dropall');
disp('TDT Initiated Successfully!');
out=1;

%% PUT A SQUARE WAVE IN DAMA BUFFER
t=(1/1000):(1/1000):200;
len=length(t);
y=square(t);
dbn=[];
cnt=1;
w=s232('freewords');
while(((cnt+1)*len)<w)
    s232('dropall')
    s232('pushf',y,len);
    s232('scale',10000);
    dbn(cnt)=s232('_allot16',len);
    s232('qpop16',dbn(cnt))
    disp(['check dama buffer ' num2str(dbn(cnt)) ' using the S232z controller''s plotting fxn']);
    cnt=cnt+1;
end
% r=rem(w,len);
% s232('dropall')
% s232('pushf',y(1:r),r);
% s232('scale',10000);
% dbn(cnt)=s232('_allot16',r);

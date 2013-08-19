function out=TDT_init
% subroutine to initiate TDT and S232 Controller

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
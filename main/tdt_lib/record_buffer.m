function [out]=record_buffer(ch, REC_BUF,DEC_BUF,TDT,signal_flag,display_flag)
% record buffers to file
% ch -- recording channel (#1,#2,etc)
% buf -- buffer indicator in record sequence (#1,#2,etc)
% outputs a trace to matlab of the last record buffer 
% if display_flag == true

S232('dropall');
S232('qpush16',REC_BUF);
S232('decimate',TDT.dec_factor);
if(display_flag>0)
	S232('qdup');
	out=S232('popf');
end
S232('make',0,signal_flag);
S232('qpop16',DEC_BUF);
S232('dama2disk16',DEC_BUF,TDT.outFN{ch},1);
S232('dropall');
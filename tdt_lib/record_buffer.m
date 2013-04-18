function [varargout]=record_buffer(REC_BUF,DEC_BUF,signal_flag,display_flag)
% record buffers to file
% ch -- recording channel (#1,#2,etc)
% buf -- buffer indicator in record sequence (#1,#2,etc)
% outputs a trace to matlab of the last record buffer 
% if display_flag == true
global TDT

S232('qpush16',REC_BUF);
S232('decimate',TDT.dec_factor);
if display_flag
	S232('qdup');
	varargout=S232('popf');
end
S232('make',0,signal_flag);
S232('qpop16',DEC_BUF);
S232('dama2disk16',DEC_BUF,TDT.outFN(ch));


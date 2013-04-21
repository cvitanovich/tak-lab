function [var_name] = dama2pc(damabuf)
%dama2pc:			Get a sound from a DAMA buffer and put into PC buffer
%damabuf:			must name an allocated DAMA integer16 buffer

if(~exist('damabuf'))
   error('DAMA buffer variable does not exist');
end

while(S232('APactive')) pause(0); end
S232('qpush16',damabuf);
while(S232('APactive')) pause(0); end
var_name = S232('pop16');   
   
return
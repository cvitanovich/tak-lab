function [playspecbuf] = make_playspecbuf(dama_bufs)
%[playspecbuf] = make_playspecbuf(dama_bufs)
%Sets up play specification for 1 or more buffers for simultaneous play

if(length(dama_bufs) < 10)
   numpush = 10;
else
   numpush = length(dama_bufs) + 2;
end


while (S232('APactive'))  pause(0); end
playspecbuf = S232('_allot16',numpush);


S232('dpush', length(dama_bufs));
S232('value',0);
for nbuf = 1:length(dama_bufs)
   S232('make',nbuf-1,dama_bufs(nbuf));
end
S232('make',length(dama_bufs),0);
S232('qpop16',playspecbuf);

return
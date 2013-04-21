function [playseqbuf] = make_playseqbuf(dama_buf)
%[playseqbuf] = make_playseqbuf(dama_buf)

while (S232('APactive'))  pause(0); end
playseqbuf = S232('_allot16',10);

S232('dpush', 10);
S232('value',0);
S232('make',0,dama_buf);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',playseqbuf);

return
function convert()
% CONVERT is a batch job to convert all result files from
%	  psychoacoustic tests to MATLAB format.

res2mat('sa_4afc.001','sa_4afc1');
res2mat('sa_4afc.002','sa_4afc2');
res2mat('sa_4afc.003','sa_4afc3');
res2mat('sb_4afc.001','sb_4afc1'); 
res2mat('sb_4afc.002','sb_4afc2'); 
res2mat('sb_4afc.003','sb_4afc3'); 
res2mat('sb_4afc.00x','sb_4afcx'); 
res2mat('sc_4afc.001','sc_4afc1'); 
res2mat('sc_4afc.002','sc_4afc2'); 
res2mat('sc_4afc.003','sc_4afc3'); 
res2mat('salocal.001','saloc01'); 
res2mat('salocal.002','saloc02'); 
res2mat('salocal.003','saloc03'); 
res2mat('sblocal.001','sbloc01'); 
res2mat('sblocal.002','sbloc02'); 
res2mat('sblocal.003','sbloc03'); 
res2mat('sclocal.001','scloc01'); 
res2mat('sclocal.002','scloc02'); 
res2mat('sclocal.003','scloc03'); 


function [] = init_tdt();

%function [] = init_tdt()
% initialize all SystemII hardware
din = 1;

if(S232('S2init', 0, 'INIT_PRIMARY', 10000) == 0)
    disp('Cannot initialize a secondary process')
    return;
else
   disp('Initialized secondary process');
end

if(S232('APlock',100, 0)==0)
   disp('FAILED to get AP lock for initialization');
   s232('S2close');
   return;
else
   disp('AP lock');
end

err = S232('getS2err');
switch err
case 0
case 1: disp('APOS error'); return;
case 2: disp('XBUS error'); return;
end


C_DATA 	= 8;
M_PULSE	= 3;
M_BIT 	= 1;

% get lock
if(S232('XBlock', 100, 0) == 0)
    disp('Cannot acquire lock on X Bus')
    s232('APunlock', 0);
    s232('S2close');
    return;
else
    disp('XBlock in Init_tdt');
end

S232('PD1clear',din);
DUALADC = s232('PD1export','DUALADC',0);

s232('PD1mode',din,DUALADC);
s232('PD1npts',din,600);
s232('PD1srate',din,200);
s232('PD1arm',din);
s232('PD1go',din);
S232('PD1stop',din);
S232('PD1clrIO',din);

return   
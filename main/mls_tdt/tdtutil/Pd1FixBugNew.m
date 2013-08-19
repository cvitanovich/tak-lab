function [] = PD1FixBugNew (din)
%function [] = PD1FixBugNew (din)

if isempty(din) din = 1;	end

C_DATA 	= 8;
M_PULSE	= 3;
M_BIT 	= 1;

% get lock
if(S232('XBlock',100, 0)==0)
   disp('No XBlock in PD1fixbugNew');
   return;
else
    disp('XBlock in PD1fixbugNew');
end


DUALADC = s232('PD1export','DUALADC',0);

s232('PD1mode',din,DUALADC);
s232('PD1npts',din,600);
s232('PD1srate',din,200);
s232('PD1arm',din);
s232('PD1go',din);

return
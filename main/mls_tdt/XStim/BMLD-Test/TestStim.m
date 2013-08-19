% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_PRIMARY', 1000) == 0)
   disp('Cannot initialize a secondary process')
   return;
end

if(S232('APlock', 100, 0) == 0)
  disp('Cannot acquire lock on AP2 Card')
  s232('S2close');
  return;
end

if(S232('XBlock', 100, 0) == 0)
  disp('Cannot acquire lock on X Bus')
  s232('APunlock', 0);
  s232('S2close');
  return;
end

%Specify DAMA buffers
stimbuf_left_1			= 1;
stimbuf_right_1		= 2;
stimbuf_left_2			= 7;
stimbuf_right_2		= 8;


stimuli_dir = 'e:\spezio\matlab\scripts\tdt\Xstim\SpaceTest\Stimuli\';
hrtfdirmat = sph2dbl(mtlrdir('e:\spezio\matlab\scripts\tdt\XStim\HRTFFilts\913f_tdAD'));
Fs = 30000;
itdfiltlen = 256;
ephonefiltlen = 255;
hrtffiltlen = 255;

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = 100; %ms
silence_len = (DUR_silence * round(Fs/1000));

%Make Stimulus buffers
DUR = 100;
space_flag = 1;
if(space_flag == 1) %ILDAlone Test
   S232('allot16',stimbuf_left_1,(silence_len + DUR*round(Fs/1000)) + itdfiltlen + ephonefiltlen + hrtffiltlen);
   S232('allot16',stimbuf_right_1,(silence_len + DUR*round(Fs/1000)) + itdfiltlen + ephonefiltlen + hrtffiltlen);
   S232('allot16',stimbuf_left_2,(silence_len + DUR*round(Fs/1000)) + itdfiltlen + ephonefiltlen + hrtffiltlen);
   S232('allot16',stimbuf_right_2,(silence_len + DUR*round(Fs/1000)) + itdfiltlen + ephonefiltlen + hrtffiltlen);
elseif(space_flag == 0) %VAS Test
   S232('allot16',stimbuf_left_1,(silence_len + DUR*round(Fs/1000)) + ephonefiltlen + hrtffiltlen);
   S232('allot16',stimbuf_right_1,(silence_len + DUR*round(Fs/1000)) + ephonefiltlen + hrtffiltlen);
   S232('allot16',stimbuf_left_2,(silence_len + DUR*round(Fs/1000)) + ephonefiltlen + hrtffiltlen);
   S232('allot16',stimbuf_right_2,(silence_len + DUR*round(Fs/1000)) + ephonefiltlen + hrtffiltlen);
end   

locscriptfid = fopen('e:\spezio\matlab\scripts\tdt\Xstim\SpaceTest\Stimuli\Checkn90.scr','r');
locs = fscanf(locscriptfid,'%i %i',[2 inf]);
for locnum = 1:length(locs)
if(hrtfdirmat(1,locnum) ~= 90)
disp(['File is stimbuf_right_' num2str(locs(1,locnum)) '_' num2str(locs(2,locnum))])   
fname = ['stimbuf_left_' num2str(locs(1,locnum)) ...
      '_' num2str(locs(2,locnum))];
evalstr = ['S232(''disk2dama16'',stimbuf_left_1,' ...
      [' ''' stimuli_dir fname ''' '] ...
      ',0);'];
eval(evalstr);
temp_left = dama2pc(stimbuf_left_1);
fname = ['stimbuf_right_' num2str(locs(1,locnum)) ...
      '_' num2str(locs(2,locnum))];
evalstr = ['S232(''disk2dama16'',stimbuf_right_1,' ...
      [' ''' stimuli_dir fname ''' '] ...
      ',0);'];
eval(evalstr);
temp_right = dama2pc(stimbuf_right_1);

plot(temp_left,'g')
hold
plot(temp_right,'r')
pause(0.2)
end
end


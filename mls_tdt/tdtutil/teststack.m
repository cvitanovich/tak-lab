clear
ephonefname = 'e:\spezio\matlab\scripts\tdt\Xstim\HRTFfilts\901etbks.imp';
Fs = 30000;

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

ephonefilt_left  = (mtlrch(ephonefname,1))';
ephonefilt_right = (mtlrch(ephonefname,2))';

t = 0:1/Fs:0.02-(1/Fs);
sig = sin(2*pi*2000*t);
sig2 = rand(1,255);


%pc2stack(sig);
%pc2stack(sig2);
pc2stack(ephonefilt_left);
temp = stack2pc;

%plot(10000*sig);
%plot(10000*sig2);
plot(10000*ephonefilt_left/max(abs(ephonefilt_left)))
hold on
plot(temp,'r');

S232('S2close') % close application and release AP2 and XBUS locks


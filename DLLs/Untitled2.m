global C_;
global M_;

run e:\kip\matlab\mls_tdt\tdtutil\globals_mii

% params
ADpts = 16000;
trigger_delay = 0;
pulse_period = 20;           % usec
dataout = zeros(ADpts,1);
ADchan = 1;

%%%%pulse_count = 260;

% init main unit
m100x( C_.INIT );

% init deep FIFO
m110dx( C_.INIT );
m110dx( C_.CLOCK, 10);
m110dx( C_.MODE, M_.PST );

input('<CR> to start clock')
%Start clock
m110dx( C_.START);
%Send pulse: PD1 GO!
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse


input('<CR> to stop clock')



m110dx( C_.STOP);
spikes = m110dx( C_.DATA, 1000); %Take 1000 spikes max

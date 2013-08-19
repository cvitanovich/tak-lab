globals_mii;

mii_us_per_sample = 10; %microsecond per sample
mii_separation = 100; %only take events separated by 100 samples (i.e., 1 ms)
for i = 1:100
    % inititialize
    m100x( C_.INIT );
    m110dx( C_.INIT );
    m110dx( C_.CLOCK, mii_us_per_sample);
    m110dx( C_.MODE, M_.PST );
   %Start clock
    m110dx( C_.START);
    %Send pulse: PD1 GO!
    m101x( C_.DATA,M_.BIT,M_.PULSE,0);                  %Use port 0 for the pulse
    pause(30);
    m110dx( C_.STOP);
         
    % read spikes
    spikes = m110dx( C_.DATA, 5000); 
    spikes = spikes(find(spikes ~= 0))                 %Get clock events that are spikes
    ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
    spikes = spikes(ind)/(1000/mii_us_per_sample);
    Nspikes = length(spikes);
    pause(1)    
end
     
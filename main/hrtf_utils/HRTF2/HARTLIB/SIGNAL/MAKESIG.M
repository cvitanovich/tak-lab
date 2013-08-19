% MAKESIG calculates signals for psychoacoustic test (pink noise, 200 ms, all directions)

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 

% 21.October 96


close all
wnoise = mkwnoise(2.0,44100);

a = [1 -0.99];

pnoise = filter(1,a,wnoise);

b1=1.100*44100+1;
e1=1.120*44100+1;
b2=1.380*44100+1;
e2=1.400*44100+1;

pnoise=fade(pnoise,b1,e1,b2,e2);
dir = 1:122;
n_dir = length(dir);

for act_dir=1:n_dir
        tfl=conv(HRIR(:,dir(act_dir)*2-1),HPIR_INV(:,1));
	left_channel = conv(pnoise,tfl);
      
        tfr=conv(HRIR(:,dir(act_dir)*2),HPIR_INV(:,2));
	right_channel = conv(pnoise,tfr);

        
        left_channel = left_channel(:);
        right_channel= right_channel(:);

        [m,n] = size(left_channel);

        outbuf = zeros(2*m,1); 
        outbuf(1:2:2*m) = left_channel;
        outbuf(2:2:2*m) = right_channel;

        m = max(outbuf);
        outbuf=outbuf/m*16000;

        outfile = ['stimulus_' int2str(act_dir)] 
        writebin(outfile,outbuf,'int16');
end;


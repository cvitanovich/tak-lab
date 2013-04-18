% ELEVSPEC  produces an image with spectral response as a function of
%           elevation for some azimuth value
%  load spectrum into freqspec1 and freqspec2
%  locations go in azim_all and elev_all

azimuth = 0;

%  filter out rows of freqspec1 that correspond with azim==0
med_spec1 = freqspec1(azim_all==0,:);
med_spec2 = freqspec2(azim_all==0,:);
med_elev = elev_all(azim_all==0);

% sort rows by elevation
[med_elev,i] = sort(med_elev);
med_spec1 = med_spec1(i',:);
med_spec2 = med_spec2(i',:);

% chop off the bottom observation
%med_elev = med_elev(1:36);
%med_spec1 = med_spec1(1:36,:);
%med_spec2 = med_spec2(1:36,:);

% now make the graphs
med_spec1 = med_spec1;
med_spec2 = med_spec2;
sf = 1000
ef = 11000
[r,c] = size(med_spec1);
freqs = (1:c)*((ef-sf)/c)+1000;

% LEFT
sh1 = surf(freqs, med_elev, med_spec1);
view(-22.5, 60);
xlabel('frequency');
ylabel('elevation');
title('left ear freq resp by elevation medial plane');

% RIGHT
figure;
sh3 = surf(freqs, med_elev, med_spec2);
view(-22.5, 60);
xlabel('frequency');
ylabel('elevation');
title('right ear freq resp by elevation medial plane');

% DIFFERENCE
diff_spec = med_spec1-med_spec2;
figure;
sh5 = surf(freqs, med_elev, diff_spec);
view(-22.5, 60);
xlabel('frequency');
ylabel('elevation');
title('binaural diff freq resp by elevation medial plane');

% pcolor plot of left
figure
%subplot(2,2,1)
sh2 = surf(freqs, med_elev, med_spec1);
colormap(jet);
view(0, 90);
set(sh2, 'EdgeColor', 'flat');
xlabel('frequency');
ylabel('elevation');
title('left freq resp by elev medial');


% pcolor plot of right
figure
%subplot(2,2,2)
sh4 = surf(freqs, med_elev, med_spec1);
view(0, 90);
set(sh4, 'EdgeColor', 'flat');
colormap(jet);
xlabel('frequency');
ylabel('elevation');
title('right freq resp by elev medial');

% pcolor plot of difference
figure
%subplot(2,2,3)
sh6 = surf(freqs, med_elev, diff_spec);
view(0, 90);
set(sh6, 'EdgeColor', 'flat');
colormap(jet);
xlabel('frequency');
ylabel('elevation');
title('diff freq resp by elev medial');

function Neuron = make_neuronstruct2(Site,...
      bird_number,...
      ILD_matrix,...
      HRTFinfo,...
      threshhold,...
      test_numbers,...
      plotflag);

%Function to create a neuron "structure" for use in further analysis
%Requires that the HRTF information already be in the workspace
%Designed to be used with new data collection method -- not David's DOS system
%The following are the properties of the structure Neuron
%  ml_position    - medial/lateral pos
%  ap_position    - anterior/posterior ps
%  depth
%  threshhold(noise)
%  sr_testpars    - Spontaneous rate test parameters
%  sr_mean        - Spontaneuos rate mean (scalar)
%  sr_std         - Spontaneous rate standard deviation (scalar)
%  abi_testpars   - Avg Binaural Intensity test parameters
%  abi_abiaxis    - ABI x-axis
%  abi_mean       - ABI mean curve
%  abi_std        - ABI standard deviation
%  itd_testpars   - Interaural Time Difference test parameters
%  itd_itdaxis    - ITD x-axis
%  itd_mean       - ITD mean curve
%  itd_std        - ITD standard deviation
%  abif_testpars  - ABI/Frequency test parameters
%  abif_freqaxis  - ABI/Frequency freq axis
%  abif_abiaxis   - ABI/Frequency ABI axis
%  abif_meansurf  - ABI/F mean response surface
%  abif_stdsurf   - ABI/F standard deviation of response surface
%  tif_testpars   - Tonal ILD/Freq test parameters
%  tif_freqaxis   - TIF frequency axis
%  tif_ildaxis    - TIF ILD axis
%  tif_meansurf   - TIF mean response surface
%  tif_stdsurf    - TIF standard deviation of response surface
%  bif_testpars   - BandPass ILD/Freq test parameters
%  bif_freqaxis   - BIF frequency axis
%  bif_ildaxis    - BIF ILD axis
%  bif_meansurf   - BIF mean response surface
%  bif_stdsurf    - BIF standard deviation of response surface

%ILDAlone - these properties each have multiple entries due to the different frequency ranges
%  ia_testpars    - ILDAlone test parameters
%  ia_locs        - ILDAlone test locations
%  ia_azi         - Azimuth for ILDAlone diamond
%  ia_ele         - Elevation for ILDAlone diamond
%  ia_meanarray   - ILDAlone mean response surface in array format
%  ia_stdarray    - ILDAlone standard deviation of response surface
%  ia_diamond     - ILDAlone diamond RS
%  tia_meanarray  - PREDICTED Tonal ILDAlone array
%  tia_azi        - Azimuth for PREDICTED Tonal ILDAlone diamond RS
%  tia_ele        - Elevation for PREDICTED Tonal ILDAlone diamond RS
%  tia_diamond    - PREDICTED Tonal ILDAlone diamond RS
%  bia_meanarray  - PREDICTED BP ILDAlone array
%  bia_azi        - Azimuth for PREDICTED BP ILDAlone diamond RS
%  bia_ele        - Elevation for PREDICTED BP ILDAlone diamond RS
%  bia_diamond    - PREDICTED BP ILDAlone diamond RS

%True Space - these properties each have multiple entries due to the different frequency ranges
%  ts_testpars    - True Space test parameters
%  ts_locs        - True Space test locations
%  ts_meanarray   - True Space mean response surface in array format
%  ts_stdarray    - True Space standard deviation of response surface
%  ts_diamond     - True Space diamond

%Specify position of test in test_numbers
SR          = 1; %Spontaneous Rate (spikes/stim)
ABL         = 2; %Average Binaural Level
ITD         = 3; %Interaural Time Difference
ABLF        = 4; %ABL-Frequency
TIF         = 5; %Tonal ILD-Frequency
TSIF			= 6; %Tone Stack ILD-Frequency
GIF         = 7; %Gammatone ILD-Frequency
GSIF			= 8; %Gammatone Stack ILD-Frequency
BBIF			= 9; %Broadband ILD-Frequency
FD				= 10; %Frequency Drop, 4 entries total
IA          = 14; %ILDAlone,  12 entries total
TS          = 26; %TrueSpace, 12 entries total

sratetest 	= test_numbers(SR);
abltest 		= test_numbers(ABL);
itdtest 		= test_numbers(ITD);
ablftest		= test_numbers(ABLF);
tiftest		= test_numbers(TIF);
tsiftest		= test_numbers(TSIF);
giftest		= test_numbers(GIF);
gsiftest		= test_numbers(GSIF);
bbiftest		= test_numbers(BBIF);
fdtest		= test_numbers(FD:FD+3);
iatest		= test_numbers(IA:IA+11);
tstest		= test_numbers(TS:TS+11);

Neuron.threshhold = threshhold;

time_window = [0 110];

%ABL data
if(abltest ~= -1)
   disp('Processing ABL test...')
   Neuron.abi_testpars = Site.test{abltest}.params;
   loabl = Site.test{abltest}.params.loabl;
   hiabl = Site.test{abltest}.params.hiabl;
   numabls = Site.test{abltest}.params.numabls;
   Neuron.abi_abiaxis = loabl:(hiabl-loabl)/(numabls-1):hiabl;
   dim2vals = [];
   reps_to_include = [];
   [Neuron.abi_mean, Neuron.abi_std] = ...
      proc_test(Site.test{abltest}.datamatrix,...
      Neuron.abi_abiaxis,...
      dim2vals,...
      time_window,...
      reps_to_include);
end

%ITD data
if(itdtest ~= -1)
   disp('Processing ITD test...')
   Neuron.itd_testpars = Site.test{itdtest}.params;
   loitd = Site.test{itdtest}.params.loitd;
   hiitd = Site.test{itdtest}.params.hiitd;
   numitds = Site.test{itdtest}.params.numitds;
   Neuron.itd_itdaxis = loitd:(hiitd-loitd)/(numitds-1):hiitd;
   dim2vals = [];
   reps_to_include = [];
   [Neuron.itd_mean, Neuron.itd_std] = ...
      proc_test(Site.test{itdtest}.datamatrix,...
      Neuron.itd_itdaxis,...
      dim2vals,...
      time_window,...
      reps_to_include);
end

%FD data
for FD_num = 1:4
   if(fdtest(FD_num) ~= -1)
      disp('Processing FD test...')
      Neuron.fd_testpars{FD_num} = Site.test{fdtest(FD_num)}.params;
      lofreq = Site.test{fdtest(FD_num)}.params.lofreq;
      hifreq = Site.test{fdtest(FD_num)}.params.hifreq;
      numfreqs = Site.test{fdtest(FD_num)}.params.numfreqs;
      Neuron.fd_freqaxis{FD_num} = lofreq:(hifreq-lofreq)/(numfreqs-1):hifreq;
      dim2vals = [];
      reps_to_include = [];
      [Neuron.fd_mean{FD_num}, Neuron.fd_std{FD_num}] = ...
         proc_test(Site.test{fdtest(FD_num)}.datamatrix,...
         Neuron.fd_freqaxis{FD_num},...
         dim2vals,...
         time_window,...
         reps_to_include);
   end
end

%ILDAlone data
for IA_num = 1:12
   if(iatest(IA_num) ~= -1)
      disp('Processing ILDAlone test...')
      Neuron.ia_testpars{IA_num} = Site.test{iatest(IA_num)}.params;
      Neuron.ia_locs{IA_num} = Site.test{iatest(IA_num)}.params.locations;
      reps_to_include = [];
      [Neuron.ia_mean{IA_num}, Neuron.ia_std{IA_num}] = ...
         proc_test_space(Site.test{iatest(IA_num)}.datamatrix,...
         Neuron.ia_locs{IA_num},...
         time_window,...
         reps_to_include);
      
      [Neuron.ia_azi{IA_num},Neuron.ia_ele{IA_num},temp] = ...
         array2diamond(Neuron.ia_mean{IA_num},Neuron.ia_locs{IA_num});
      
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(Neuron.ia_azi{IA_num}, Neuron.ia_ele{IA_num});
      
      % generate mask for missing points
      missmask = NaN*ones(size(temp));
      
      i = 1;
      for az = -90:5:90;
         for el = -90+abs(az):5:90-abs(az)
            if (~(Neuron.ia_locs{IA_num}(1,:)==el & Neuron.ia_locs{IA_num}(2,:)==az))
               missmask(AZ==az & EL==el) = 1;
            else
               missmask(AZ==az & EL==el) = 0;
            end; 
         end;
      end;
      
      % replace all valid locations presently containing NaN's with zeros
      ri = (missmask==1);
      temp(ri) = zeros(size(temp(ri)));
      % generate interpolation function - uses surrounding 4 squares to generate missing square
      intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
      temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
      intval_old = conv2(temp, intfun,'same')/4;
      intval = intval_old.*missmask;
      Neuron.ia_diamond{IA_num} = intval + temp;
   end
end

%True Space data
for TS_num = 1:12
   if(tstest(TS_num) ~= -1)
      disp('Processing True Space test...')
      Neuron.ts_testpars{TS_num} = Site.test{tstest(TS_num)}.params;
      Neuron.ts_locs{TS_num} = Site.test{tstest(TS_num)}.params.locations;
      reps_to_include = [];
      [Neuron.ts_mean{TS_num}, Neuron.ts_std{TS_num}] = ...
         proc_test_space(Site.test{tstest(TS_num)}.datamatrix,...
         Neuron.ts_locs{TS_num},...
         time_window,...
         reps_to_include);
      
      [Neuron.ts_azi{TS_num},Neuron.ts_ele{TS_num},temp] = ...
         array2diamond(Neuron.ts_mean{TS_num},Neuron.ts_locs{TS_num});
      
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(Neuron.ts_azi{TS_num}, Neuron.ts_ele{TS_num});
      
      % generate mask for missing points
      missmask = NaN*ones(size(temp));
      
      i = 1;
      for az = -90:5:90;
         for el = -90+abs(az):5:90-abs(az)
            if (~(Neuron.ts_locs{TS_num}(1,:)==el & Neuron.ts_locs{TS_num}(2,:)==az))
               missmask(AZ==az & EL==el) = 1;
            else
               missmask(AZ==az & EL==el) = 0;
            end; 
         end;
      end;
      
      % replace all valid locations presently containing NaN's with zeros
      ri = (missmask==1);
      temp(ri) = zeros(size(temp(ri)));
      % generate interpolation function - uses surrounding 4 squares to generate missing square
      intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
      temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
      intval_old = conv2(temp, intfun,'same')/4;
      intval = intval_old.*missmask;
      Neuron.ts_diamond{TS_num} = intval + temp;
   end
end

disp(' ')

if(iatest(1) ~= -1)
   %Get those locations in the ILDmatrix that match the measured ILDAlone locations
   for num_loc = 1:size(Neuron.ia_locs{1},2)
      ILDmat_index(num_loc) =...
         max(find(HRTFinfo.location_matrix(1,:) == Neuron.ia_locs{1}(1,num_loc) &...
         HRTFinfo.location_matrix(2,:) == Neuron.ia_locs{1}(2,num_loc)));
   end
   ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   %Find the best ILDAlone RS location
   IA_bestloc_ind = max(find(Neuron.ia_mean{1} == ...
      max(Neuron.ia_mean{1})));
   [y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - Site.test{iatest(1)}.params.lofreq));
   [y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - Site.test{iatest(1)}.params.hifreq));
end

if(iatest(1) ~= -1 & tstest(1) ~= -1)
   %Find the best True Space RS location
   TS_bestloc_ind = max(find(Neuron.ts_mean{1} == ...
      max(Neuron.ts_mean{1})));
end

%TIF data
if(tiftest ~= -1)
   disp('Processing TIF test...')
   Neuron.tif_testpars = Site.test{tiftest}.params;
   lofreq = Site.test{tiftest}.params.lofreq;
   hifreq = Site.test{tiftest}.params.hifreq;
   numfreqs = Site.test{tiftest}.params.numfreqs;
   Neuron.tif_freqaxis = lofreq:(hifreq-lofreq)/(numfreqs-1):hifreq;
   loild = Site.test{tiftest}.params.loild;
   hiild = Site.test{tiftest}.params.hiild;
   numilds = Site.test{tiftest}.params.numilds;
   Neuron.tif_ildaxis = loild:(hiild-loild)/(numilds-1):hiild;
   reps_to_include = [];
   [Neuron.tif_mean, Neuron.tif_std] = ...
      proc_test(Site.test{tiftest}.datamatrix,...
      Neuron.tif_freqaxis,...
      Neuron.tif_ildaxis,...
      time_window,...
      reps_to_include);
   
   if(Site.test{tiftest}.params.loc_flag == 1)
      disp('Correcting for optimal spatial filter...')
      %Get info for optimal spatial filter
      azel = Site.test{tiftest}.params.loc_azel;
      optlocnum = find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
         HRTFinfo.location_matrix(2,:) == azel(1));
      [Neuron] = iftransform_5(Neuron,ILD_matrix(:,optlocnum),HRTFinfo,'tif');
   end
   
   disp('Transforming into VAS...')
   [Neuron.tia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron.tif_mean,...
      Neuron.tif_freqaxis,...
      Neuron.tif_ildaxis);
   [Neuron.tia_azi,Neuron.tia_ele,Neuron.tia_diamond] = ...
      array2diamond(Neuron.tia_meanarray,HRTFinfo.location_matrix);
end

%GIF data
if(giftest ~= -1)
   disp('Processing GIF test...')
   Neuron.gif_testpars = Site.test{giftest}.params;
   lofreq = Site.test{giftest}.params.lofreq;
   hifreq = Site.test{giftest}.params.hifreq;
   numfreqs = Site.test{giftest}.params.numfreqs;
   Neuron.gif_freqaxis = lofreq:(hifreq-lofreq)/(numfreqs-1):hifreq;
   loild = Site.test{giftest}.params.loild;
   hiild = Site.test{giftest}.params.hiild;
   numilds = Site.test{giftest}.params.numilds;
   Neuron.gif_ildaxis = loild:(hiild-loild)/(numilds-1):hiild;
   reps_to_include = [];
   [Neuron.gif_mean, Neuron.gif_std] = ...
      proc_test(Site.test{giftest}.datamatrix,...
      Neuron.gif_freqaxis,...
      Neuron.gif_ildaxis,...
      time_window,...
      reps_to_include);
   
   if(Site.test{giftest}.params.loc_flag == 1)
      disp('Correcting for optimal spatial filter...')
      %Get info for optimal spatial filter
      azel = Site.test{giftest}.params.loc_azel;
      optlocnum = max(find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
         HRTFinfo.location_matrix(2,:) == azel(1)));
      [Neuron] = iftransform_5(Neuron,ILD_matrix(:,optlocnum),HRTFinfo,'gif');
   end
   
   disp('Transforming into VAS...')
   [Neuron.gia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron.gif_mean,...
      Neuron.gif_freqaxis,...
      Neuron.gif_ildaxis);
   [Neuron.gia_azi,Neuron.gia_ele,Neuron.gia_diamond] = ...
      array2diamond(Neuron.gia_meanarray,HRTFinfo.location_matrix);
end

%TSIF data
if(tsiftest ~= -1)
   disp('Processing TSIF test...')
   Neuron.tsif_testpars = Site.test{tsiftest}.params;
   lofreq = Site.test{tsiftest}.params.lofreq;
   hifreq = Site.test{tsiftest}.params.hifreq;
   numfreqs = Site.test{tsiftest}.params.numfreqs;
   Neuron.tsif_freqaxis = lofreq:(hifreq-lofreq)/(numfreqs-1):hifreq;
   loild = Site.test{tsiftest}.params.loild;
   hiild = Site.test{tsiftest}.params.hiild;
   numilds = Site.test{tsiftest}.params.numilds;
   Neuron.tsif_ildaxis = loild:(hiild-loild)/(numilds-1):hiild;
   reps_to_include = [];
   [Neuron.tsif_mean, Neuron.tsif_std] = ...
      proc_test(Site.test{tsiftest}.datamatrix,...
      Neuron.tsif_freqaxis,...
      Neuron.tsif_ildaxis,...
      time_window,...
      reps_to_include);
   
   if(Site.test{tsiftest}.params.loc_flag == 1)
      disp('Correcting for optimal spatial filter...')
      %Get info for optimal spatial filter
      azel = Site.test{tsiftest}.params.loc_azel;
      optlocnum = max(find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
         HRTFinfo.location_matrix(2,:) == azel(1)));
      [Neuron] = iftransform_5(Neuron,ILD_matrix(:,optlocnum),HRTFinfo,'tsif');
   end
   
   disp('Transforming into VAS...')
   [Neuron.tsia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron.tsif_mean,...
      Neuron.tsif_freqaxis,...
      Neuron.tsif_ildaxis);
   [Neuron.tsia_azi,Neuron.tsia_ele,Neuron.tsia_diamond] = ...
      array2diamond(Neuron.tsia_meanarray,HRTFinfo.location_matrix);
end

%GSIF data
if(gsiftest ~= -1)
   disp('Processing GSIF test...')
   Neuron.gsif_testpars = Site.test{gsiftest}.params;
   lofreq = Site.test{gsiftest}.params.lofreq;
   hifreq = Site.test{gsiftest}.params.hifreq;
   numfreqs = Site.test{gsiftest}.params.numfreqs;
   Neuron.gsif_freqaxis = lofreq:(hifreq-lofreq)/(numfreqs-1):hifreq;
   loild = Site.test{gsiftest}.params.loild;
   hiild = Site.test{gsiftest}.params.hiild;
   numilds = Site.test{gsiftest}.params.numilds;
   Neuron.gsif_ildaxis = loild:(hiild-loild)/(numilds-1):hiild;
   reps_to_include = [];
   [Neuron.gsif_mean, Neuron.gsif_std] = ...
      proc_test(Site.test{gsiftest}.datamatrix,...
      Neuron.gsif_freqaxis,...
      Neuron.gsif_ildaxis,...
      time_window,...
      reps_to_include);
   
   %azel = Site.test{gsiftest}.params.loc_azel;
   %optlocnum = find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
   %   HRTFinfo.location_matrix(2,:) == azel(1));
   
   if(Site.test{gsiftest}.params.loc_flag == 1)
      disp('Correcting for optimal spatial filter...')
      %Get info for optimal spatial filter
      azel = Site.test{gsiftest}.params.loc_azel;
      optlocnum = max(find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
         HRTFinfo.location_matrix(2,:) == azel(1)));
      [Neuron] = iftransform_5(Neuron,ILD_matrix(:,optlocnum),HRTFinfo,'gsif',...
         ILD_matrix_focus(:,IA_bestloc_ind));
   end
   
   disp('Transforming into VAS...')
   [Neuron.gsia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron.gsif_mean,...
      Neuron.gsif_freqaxis,...
      Neuron.gsif_ildaxis);
   [Neuron.gsia_azi,Neuron.gsia_ele,Neuron.gsia_diamond] = ...
      array2diamond(Neuron.gsia_meanarray,HRTFinfo.location_matrix);
end

%BBIF
if(bbiftest ~= -1)
   disp('Processing BBIF test...')
   Neuron.bbif_testpars = Site.test{bbiftest}.params;
   lofreq = Site.test{bbiftest}.params.lofreq;
   hifreq = Site.test{bbiftest}.params.hifreq;
   numfreqs = Site.test{bbiftest}.params.numfreqs;
   Neuron.bbif_freqaxis = lofreq:(hifreq-lofreq)/(numfreqs-1):hifreq;
   loild = Site.test{bbiftest}.params.loild;
   hiild = Site.test{bbiftest}.params.hiild;
   numilds = Site.test{bbiftest}.params.numilds;
   Neuron.bbif_ildaxis = loild:(hiild-loild)/(numilds-1):hiild;
   reps_to_include = [];
   [Neuron.bbif_mean, Neuron.bbif_std] = ...
      proc_test(Site.test{bbiftest}.datamatrix,...
      Neuron.bbif_freqaxis,...
      Neuron.bbif_ildaxis,...
      time_window,...
      reps_to_include);
   
   %azel = Site.test{bbiftest}.params.loc_azel;
   %optlocnum = find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
   %   HRTFinfo.location_matrix(2,:) == azel(1));
   
   if(Site.test{bbiftest}.params.loc_flag == 1)
      disp('Correcting for optimal spatial filter...')
      %Get info for optimal spatial filter
      azel = Site.test{bbiftest}.params.loc_azel;
      optlocnum = max(find(HRTFinfo.location_matrix(1,:) == azel(2) & ...
         HRTFinfo.location_matrix(2,:) == azel(1)));
      [Neuron] = iftransform_5(Neuron,ILD_matrix(:,optlocnum),HRTFinfo,'bbif',...
         ILD_matrix_focus(:,IA_bestloc_ind));
   end
   
   disp('Transforming into VAS...')
   [Neuron.bbia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron.bbif_mean,...
      Neuron.bbif_freqaxis,...
      Neuron.bbif_ildaxis);
   [Neuron.bbia_azi,Neuron.bbia_ele,Neuron.bbia_diamond] = ...
      array2diamond(Neuron.bbia_meanarray,HRTFinfo.location_matrix);
end


if(plotflag)
   hfig = figure;
   set(hfig,'Units','inches',...
      'Position',[0 0 7.5 10]);
   
   %List information for the Neuron
   Neuron_info = [num2str(Site.test{iatest(1)}.params.bird_number) ...
         Site.test{iatest(1)}.params.side_of_brain ...
         num2str(Site.test{iatest(1)}.params.session_num) ','...
         ' Rec Site '...
         num2str(Site.test{iatest(1)}.params.recording_site) ','...
         ' Depth = '...
         num2str(Site.test{iatest(1)}.params.depth) ','...
         ' Stim Dur = '...
         num2str(Site.test{iatest(1)}.params.curr_stimdur) ','...
         ' ISI = '...
         num2str(Site.test{iatest(1)}.params.test_ISI)];
   htext = uicontrol('Parent',hfig,...
      'Style','text',...
      'Units','normal',...
      'Position',[0 0.98 1 0.02],...
      'String',Neuron_info,...
      'FontSize',10,...
      'FontWeight','bold');
   
   %Set horizontal figure position
   horiz_pos =		[0.75 0.50 0.25 0.05] ;
   subplot_ht = 	[0.15 0.15 0.20 0.15];
   
   %GIF
   if(giftest ~= -1)
      hgif = subplot('Position',[0.05 horiz_pos(1) 0.4 subplot_ht(1)]);
      set(hgif,'FontSize',8);
      plotsurf(Neuron.gif_freqaxis,Neuron.gif_ildaxis,Neuron.gif_mean'/max(max(Neuron.gif_mean)));
      if(iatest(1) ~= -1)
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
            'm',...
            'LineWidth',1.5);
      end
      if(tstest(1) ~= -1)
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
            'w',...
            'LineWidth',1.5);
      end
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      title('ILD-Freq, Gammatones');
   end
   
   %TSIF
   if(tsiftest ~= -1)
      htsif = subplot('Position',[0.55 horiz_pos(1) 0.4 subplot_ht(1)]);
      set(htsif,'FontSize',8);
      plotsurf(Neuron.tsif_freqaxis,Neuron.tsif_ildaxis,Neuron.tsif_mean'/max(max(abs(Neuron.tsif_mean))));
      hold on
      plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
         ILD_matrix(minfreqind:maxfreqind,optlocnum),...
         'r',...
         'Linewidth',2);
      if(iatest(1) ~= -1)
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
            'm-.',...
            'LineWidth',1.5);
      end
      if(tstest(1) ~= -1)
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
            'w--',...
            'LineWidth',1.5);
      end
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      title('ILD-Freq, Tone Stack');
      
      
      if(isfield(Neuron,'tsif_Zscores'))
         hztsif = subplot('Position',[0.05 horiz_pos(2) 0.4 subplot_ht(2)]);
         set(hztsif,'FontSize',8);
         plotsurf(Neuron.tsif_freqaxis,Neuron.tsif_ildaxis,Neuron.tsif_Zscores');
         xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
         title('Z-scores for Tone Stack Data');
         colorbar
      end
   end
   
   %GSIF
   if(bbiftest == -1 & gsiftest ~= -1)
      hgsif = subplot('Position',[0.55 horiz_pos(1) 0.4 subplot_ht(1)]);
      set(hgsif,'FontSize',8);
      plotsurf(Neuron.gsif_freqaxis,Neuron.gsif_ildaxis,Neuron.gsif_mean'/max(max(abs(Neuron.gsif_mean))));
      hold on
      plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
         ILD_matrix(minfreqind:maxfreqind,optlocnum),...
         'r',...
         'Linewidth',2);
      if(iatest(1) ~= -1)
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
            'm-.',...
            'LineWidth',1.5);
      end
      if(tstest(1) ~= -1)
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
            'w--',...
            'LineWidth',1.5);
      end
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      title('ILD-Freq, Gammatone Stack');
      
      if(isfield(Neuron,'gsif_orig_meansurf'))
         hbbiforig = subplot('Position',[0.55 horiz_pos(2) 0.4 subplot_ht(2)]);
         set(hbbiforig,'FontSize',8);
         plotsurf(Neuron.gsif_freqaxis,Neuron.gsif_ildaxis,Neuron.gsif_orig_meansurf'/max(max(abs(Neuron.gsif_orig_meansurf))));
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix(minfreqind:maxfreqind,optlocnum),...
            'r',...
            'Linewidth',2);
         if(iatest(1) ~= -1)
            hold on
            plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
               ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
               'm-.',...
               'LineWidth',1.5);
         end
         if(tstest(1) ~= -1)
            plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
               ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
               'w--',...
               'LineWidth',1.5);
         end
         colorbar
         xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
         title('Pre-Transformed ILD-Freq, BBNoise');
      end
      
      if(isfield(Neuron,'gsif_Zscores'))
         hzgsif = subplot('Position',[0.05 horiz_pos(2) 0.4 subplot_ht(2)]);
         plotsurf(Neuron.gsif_freqaxis,Neuron.gsif_ildaxis,Neuron.gsif_Zscores');
         set(hzgsif,'FontSize',8);
         xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
         title('Z-scores for Gammatone Stack Data');
         colorbar
      end
   end
   
   %BBIF
   if(bbiftest ~= -1)
      hbbif = subplot('Position',[0.55 horiz_pos(1) 0.4 subplot_ht(1)]);
      set(hbbif,'FontSize',8);
      plotsurf(Neuron.bbif_freqaxis,Neuron.bbif_ildaxis,Neuron.bbif_mean'/max(max(abs(Neuron.bbif_mean))));
      hold on
      plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
         ILD_matrix(minfreqind:maxfreqind,optlocnum),...
         'r',...
         'Linewidth',2);
      if(iatest(1) ~= -1)
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
            'm-.',...
            'LineWidth',1.5);
      end
      if(tstest(1) ~= -1)
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
            'w--',...
            'LineWidth',1.5);
      end
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      title('ILD-Freq, BBNoise');
      
      if(isfield(Neuron,'bbif_orig_meansurf'))
         hbbiforig = subplot('Position',[0.55 horiz_pos(2) 0.4 subplot_ht(2)]);
         set(hbbiforig,'FontSize',8);
         plotsurf(Neuron.bbif_freqaxis,Neuron.bbif_ildaxis,Neuron.bbif_orig_meansurf'/max(max(abs(Neuron.bbif_orig_meansurf))));
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix(minfreqind:maxfreqind,optlocnum),...
            'r',...
            'Linewidth',2);
         if(iatest(1) ~= -1)
            hold on
            plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
               ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
               'm-.',...
               'LineWidth',1.5);
         end
         if(tstest(1) ~= -1)
            plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
               ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
               'w--',...
               'LineWidth',1.5);
         end
         colorbar
         xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
         title('Pre-Transformed ILD-Freq, BBNoise');
      end
      
      if(isfield(Neuron,'bbif_Zscores'))
         hzbbif = subplot('Position',[0.05 horiz_pos(2) 0.4 subplot_ht(2)]);
         plotsurf(Neuron.bbif_freqaxis,Neuron.bbif_ildaxis,Neuron.bbif_Zscores');
         set(hzbbif,'FontSize',8);
         xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
         title('Z-scores for BBNoise Data');
         colorbar
      end
   end
   
   
   %TS
   if(tstest(1) ~= -1)
      hts = subplot('Position',[0.03 horiz_pos(3) 0.2 subplot_ht(3)]);
      set(hts,'FontSize',8);
      plotdiam(Neuron.ts_azi{1},Neuron.ts_ele{1},Neuron.ts_diamond{1});
      hc = colorbar;
      set(hc,'FontSize',8);
      title('True Space');
   end
   
   %IA
   if(iatest(1) ~= -1)
      hia = subplot('Position',[0.28 horiz_pos(3) 0.2 subplot_ht(3)]);
      set(hia,'FontSize',8);
      plotdiam(Neuron.ia_azi{1},Neuron.ia_ele{1},Neuron.ia_diamond{1});
      hc = colorbar
      set(hc,'FontSize',8);
      title('ILDAlone');
   end
   
   %GIA
   if(giftest ~= -1)
      hgia = subplot('Position',[0.53 horiz_pos(3) 0.2 subplot_ht(3)]);
      set(hgia,'FontSize',8);
      plotdiam(Neuron.gia_azi,Neuron.gia_ele,Neuron.gia_diamond);
      hc = colorbar;
      set(hc,'FontSize',8);
      title('IA: Gammatones');
   end
   
   if(giftest ~= -1)
      hgia_corr = subplot('Position',[0.53 horiz_pos(4) 0.15 subplot_ht(4)]);
      set(hgia_corr,'FontSize',8);
      temp_var = (corrcoef(Neuron.ia_mean{1},Neuron.gia_meanarray(ILDmat_index))).^2;
      plot(Neuron.ia_mean{1},Neuron.gia_meanarray(ILDmat_index),'.');
      axis square
      title(str2mat('Frac Var, Gammatone IA',...
         ['r^{2} = ' num2str(temp_var(1,2))]));
   end
   
   %TSIA
   if(tsiftest ~= -1)
      htsia = subplot('Position',[0.78 horiz_pos(3) 0.2 subplot_ht(3)]);
      set(htsia,'FontSize',8);
      plotdiam(Neuron.tsia_azi,Neuron.tsia_ele,Neuron.tsia_diamond);
      hc = colorbar;
      set(hc,'FontSize',8);
      title('IA: Tone Stack');
   end
   
   if(tsiftest ~= -1)
      htsia_corr = subplot('Position',[0.78 horiz_pos(4) 0.15 subplot_ht(4)]);
      set(htsia_corr,'FontSize',8);
      temp_var = (corrcoef(Neuron.ia_mean{1},Neuron.tsia_meanarray(ILDmat_index))).^2;
      plot(Neuron.ia_mean{1},Neuron.tsia_meanarray(ILDmat_index),'.');
      axis square
      title(str2mat('Frac Var, Tone Stack IA',...
         ['r^{2} = ' num2str(temp_var(1,2))]));
   end
   
   %GSIA
   if(bbiftest == -1 & gsiftest ~= -1)
      hgsia = subplot('Position',[0.78 horiz_pos(3) 0.2 subplot_ht(4)]);
      set(hgsia,'FontSize',8);
      plotdiam(Neuron.gsia_azi,Neuron.gsia_ele,Neuron.gsia_diamond);
      hc = colorbar;
      set(hc,'FontSize',8);
      title('IA: Gammatone Stack');
   end
   
   if(tsiftest == -1 & gsiftest ~= -1)
      hgsia_corr = subplot('Position',[0.78 horiz_pos(4) 0.15 subplot_ht(4)]);
      set(hgsia_corr,'FontSize',8);
      temp_var = (corrcoef(Neuron.ia_mean{1},Neuron.gsia_meanarray(ILDmat_index))).^2;
      plot(Neuron.ia_mean{1},Neuron.gsia_meanarray(ILDmat_index),'.');
      axis square
      title(str2mat('Frac Var, Gammatone Stack IA',...
         ['      r^{2} = ' num2str(temp_var(1,2))]));
   end
   
   %BBIA
   if(bbiftest ~= -1)
      hbbia = subplot('Position',[0.78 horiz_pos(3) 0.2 subplot_ht(3)]);
      set(hbbia,'FontSize',8);
      plotdiam(Neuron.bbia_azi,Neuron.bbia_ele,Neuron.bbia_diamond);
      hc = colorbar;
      set(hc,'FontSize',8);
      title('IA: BBNoise');
   end
   
   if(bbiftest ~= -1)
      hbbia_corr = subplot('Position',[0.78 horiz_pos(4) 0.15 subplot_ht(4)]);
      set(hbbia_corr,'FontSize',8);
      temp_var = (corrcoef(Neuron.ia_mean{1},Neuron.bbia_meanarray(ILDmat_index))).^2;
      plot(Neuron.ia_mean{1},Neuron.bbia_meanarray(ILDmat_index),'.');
      axis square
      title(str2mat('Frac Var, BBNoise IA',...
         ['      r^{2} = ' num2str(temp_var(1,2))]));
   end
   
end
      
   
   
return;



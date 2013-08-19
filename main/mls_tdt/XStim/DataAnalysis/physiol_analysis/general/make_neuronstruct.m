function Neuron = make_neuronstruct(bird_number,side_of_brain,test_numbers,ILD_matrix,HRTFinfo,reps)
%Function to create a neuron "structure" for use in further analysis
%Requires that the HRTF information already be in the workspace
%and that raw processing is already done
%Test_Numbers format:
%[Threshhold SR ABI ITD ABI/F TIF BIF IA(different frequency ranges, 10 entries) TS(different ranges, 10 entries)]
%"reps" is a vector with the number of reps for each test to include
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

if(nargin < 6)
   reps = zeros(1,length(test_numbers));
end

%Specify position of test in test_numbers
Threshhold  = 1;
SR          = 2;
ABI         = 3;
ITD         = 4;
ABIF        = 5;
TIF         = 6;
BIF         = 7;
IA          = 8;  %12 entries total
TS          = 20; %12 entries total

Neuron.threshhold = test_numbers(Threshhold);

%Spontaneous rate data
if(test_numbers(SR) ~= -1)
   repinc = dezero(reps(:,SR));
   [mean,...
         std,...
         x,...
         axis,...
         Neuron.sr_testpars, ...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(SR), 1, 0, [], [], repinc);
   Neuron.sr_mean = mean(mean);
   Neuron.sr_std  = mean(std);
end

%ABI data
if(test_numbers(ABI) ~= -1)
   repinc = dezero(reps(:,ABI));
   [Neuron.abi_mean,...
         Neuron.abi_std,...
         x,...
         Neuron.abi_abiaxis,...
         Neuron.abi_testpars, ...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(ABI), 1, 0, [], [], repinc);
end

%ITD data
if(test_numbers(ITD) ~= -1)
   repinc = dezero(reps(:,ITD));
   [Neuron.itd_mean,...
         Neuron.itd_std,...
         x,...
         Neuron.itd_itdaxis,...
         Neuron.itd_testpars, ...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(ITD), 1, 0, [], [], repinc);
end

%ABIF data
if(test_numbers(ABIF) ~= -1)
   repinc = dezero(reps(:,ABIF));
   [Neuron.abif_mean,...
         Neuron.abif_std,...
         Neuron.abif_abiaxis,...
         Neuron.abif_freqaxis,...
         Neuron.abif_testpars, ...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(ABIF), 1, 0, [], [], repinc);
end

%TIF data
if(test_numbers(TIF) ~= -1)
   repinc = dezero(reps(:,TIF));
   [Neuron.tif_meansurf,...
         Neuron.tif_stdsurf,...
         Neuron.tif_ildaxis,...
         Neuron.tif_freqaxis,...
         Neuron.tif_testpars,...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(TIF), 1, 0, [], [], repinc);
   [Neuron.tia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs, HRTFinfo.location_matrix, Neuron.tif_meansurf, Neuron.tif_freqaxis, Neuron.tif_ildaxis);
   [Neuron.tia_azi,Neuron.tia_ele,Neuron.tia_diamond] = ...
      array2diamond(Neuron.tia_meanarray,HRTFinfo.location_matrix);
end

%BIF data
if(test_numbers(BIF) ~= -1)
   repinc = dezero(reps(:,BIF));
   [Neuron.bif_meansurf,...
         Neuron.bif_stdsurf,...
         Neuron.bif_ildaxis,...
         Neuron.bif_freqaxis,...
         Neuron.bif_testpars,...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(BIF), 1, 0, [], [], repinc);
   [Neuron.bia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs, HRTFinfo.location_matrix, Neuron.bif_meansurf, Neuron.bif_freqaxis, Neuron.bif_ildaxis);
   [Neuron.bia_azi,Neuron.bia_ele,Neuron.bia_diamond] = ...
      array2diamond(Neuron.bia_meanarray,HRTFinfo.location_matrix);
end

%IA data
for IA_num = 1:12
   if(test_numbers((IA-1) + IA_num) ~= -1)
      repinc = dezero(reps(:,(IA-1) + IA_num));
      %Get the mean response surface of the ILD Alone data
      [Neuron.ia_meansurf{IA_num},...
            Neuron.ia_stdsurf{IA_num},...
            dim2vals,...
            dim1vals,...
            Neuron.ia_testpars{IA_num},...
            spont_spikes, spont_dur, nincl_reps,...
            Neuron.ia_locs{IA_num},Neuron.ia_header{IA_num}] = ...
      proc_test899(bird_number, side_of_brain, test_numbers((IA-1) + IA_num),1,0,[],[],repinc);
      [Neuron.ia_azi{IA_num},Neuron.ia_ele{IA_num},temp] = ...
         array2diamond(Neuron.ia_meansurf{IA_num},Neuron.ia_locs{IA_num}');
      %Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Neuron.ia_azi{IA_num}, Neuron.ia_ele{IA_num});

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(Neuron.ia_locs{IA_num}(:,1)==el & Neuron.ia_locs{IA_num}(:,2)==az))
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
end %end if statement

end %end loop over IA experiments

%TS data
for TS_num = 1:12
   if(test_numbers((TS-1) + TS_num) ~= -1)
      repinc = dezero(reps(:,(IA-1) + IA_num));
      %Get the mean response surface of the ILD Alone data
      [Neuron.ts_meansurf{TS_num},...
            Neuron.ts_stdsurf{TS_num},...
            dim2vals,...
            dim1vals,...
            Neuron.ts_testpars{TS_num},...
            spont_spikes, spont_dur, nincl_reps,...
            Neuron.ts_locs{TS_num},Neuron.ts_header{TS_num}] = ...
      proc_test899(bird_number, side_of_brain, test_numbers((TS-1) + TS_num),1,0,[],[],repinc);
      [Neuron.ts_azi{TS_num},Neuron.ts_ele{TS_num},temp] = ...
         array2diamond(Neuron.ts_meansurf{TS_num},Neuron.ts_locs{TS_num}');
      %Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Neuron.ts_azi{TS_num}, Neuron.ts_ele{TS_num});

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(Neuron.ts_locs{TS_num}(:,1)==el & Neuron.ts_locs{TS_num}(:,2)==az))
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
end %end if statement

end %end loop over TS experiments




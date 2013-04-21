clear
%close all
%Load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ILD_matrix ITD_matrix ABI_matrix HRTFinfo test_numbers reps
cmap = 1-gray;

test_numbers = [
%[Thr  SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]
  -57  -1  -1  65    -1  66  68  70  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  71  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -62  -1  74  75    -1  76  77  79  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  80  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1  -1  85    -1  87  88  89  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  90  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1  -1  91    -1  92  93  94  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  95  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -80  -1  96  97    -1  98 100 101  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1 103 104    -1 105 106 107  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 108 109    -1 110 111 112  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 113  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 121 122    -1 123 124 125  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 126  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1  -1 127    -1 129 130 131  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 132  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 133 134    -1 135 136 137  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 138  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 141 142    -1 143 144 145  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 146  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 152 153    -1 162 158 159  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 161  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 165 166    -1 167 168 169  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 170  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 171 172    -1 173 174 175  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 176  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 206 208    -1 209 210 211  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 207 212 213  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70 220 214 216   223 217 218 219  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 221  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 228 230   234 232 233 231  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 229  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 242 248    -1 249 250 247  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 246  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 264 265    -1 274  -1 267 268  -1 269  -1 270 271 272 273   -1   -1   -1 266  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 289 290    -1 300  -1 293 294  -1 295  -1 296 297 298 299   -1   -1   -1 292  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -60  -1 318 319    -1  -1  -1 321 322  -1 323  -1  -1  -1  -1  -1   -1   -1   -1 320  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 326 327    -1 332  -1 329 330  -1 331 333 334 335  -1  -1   -1   -1   -1 328  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 337 338    -1 343  -1 340 341  -1 342  -1  -1  -1  -1  -1   -1   -1   -1 339  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 347 348    -1 356  -1 351 352  -1 353  -1  -1 355 358 357   -1   -1   -1 349  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
];

%Parameters
bird_number = 899;
side_of_brain = 'r';

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;

%Get the HRTF data
[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = ...
   get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);


%Set the time limits
time_limits = [108 138
   138 168
   168 198
   198 228
	228 258]; %ms after stim onset

for cell_num = 1:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   fig_ht = 9.5; fig_wd = 7.5;      
   mainfig = figure; set(mainfig,'Units','inches'); 
   set(mainfig,'Position',[0 0 fig_wd fig_ht]);

   for time_limit = 1:size(time_limits,1)
   
   Neuron_timelim1{cell_num}.time_limits = time_limits(time_limit,:);
   
      %Find the best ILDAlone RS location
      IA_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ia_meansurf{1} == ...
         max(Neuron{cell_num}.ia_meansurf{1})));
      if(cell_num == 1)  IA_bestloc_ind(cell_num) = 129; end
      if(cell_num == 4)  IA_bestloc_ind(cell_num) = 198; end
      if(cell_num == 7)  IA_bestloc_ind(cell_num) = 163; end
      if(cell_num == 9)  IA_bestloc_ind(cell_num) = 282; end
      if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
      if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
      if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
      if(cell_num == 18) IA_bestloc_ind(cell_num) =  78; end
      if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
      
      
		%Get those locations in the ILDmatrix that match the measured ILDAlone locations
      for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
         ILDmat_index(num_loc) =...
            max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
            HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
      end
      ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   
   %Get the tests under specified time limits
   
   %1. Tonal
   
   [Neuron_timelim1{cell_num}.tif_meansurf,...
         Neuron_timelim1{cell_num}.tif_stdsurf,...
         Neuron_timelim1{cell_num}.tif_ildaxis,...
         Neuron_timelim1{cell_num}.tif_freqaxis,...
         Neuron_timelim1{cell_num}.tif_testpars,...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number,...
      side_of_brain,...
      test_numbers(cell_num,6),...
      1,...
      0,...
      time_limits(time_limit,1),...
      time_limits(time_limit,2));
   
   [Neuron_timelim1{cell_num}.tia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron_timelim1{cell_num}.tif_meansurf,...
      Neuron_timelim1{cell_num}.tif_freqaxis,...
      Neuron_timelim1{cell_num}.tif_ildaxis);
   
   [Neuron_timelim1{cell_num}.tia_azi,...
         Neuron_timelim1{cell_num}.tia_ele,...
         Neuron_timelim1{cell_num}.tia_diamond] = ...
      array2diamond(Neuron_timelim1{cell_num}.tia_meanarray,HRTFinfo.location_matrix);
   
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
   %2. BP
   
   [Neuron_timelim1{cell_num}.bif_meansurf,...
         Neuron_timelim1{cell_num}.bif_stdsurf,...
         Neuron_timelim1{cell_num}.bif_ildaxis,...
         Neuron_timelim1{cell_num}.bif_freqaxis,...
         Neuron_timelim1{cell_num}.bif_testpars,...
         spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number,...
      side_of_brain,...
      test_numbers(cell_num,7),...
      1,...
      0,...
      time_limits(time_limit,1),...
      time_limits(time_limit,2));
   
   [Neuron_timelim1{cell_num}.bia_meanarray,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      HRTFinfo.hrtf_freqs,...
      HRTFinfo.location_matrix,...
      Neuron_timelim1{cell_num}.bif_meansurf,...
      Neuron_timelim1{cell_num}.bif_freqaxis,...
      Neuron_timelim1{cell_num}.bif_ildaxis);
   
   [Neuron_timelim1{cell_num}.bia_azi,...
         Neuron_timelim1{cell_num}.bia_ele,...
         Neuron_timelim1{cell_num}.bia_diamond] = ...
      array2diamond(Neuron_timelim1{cell_num}.bia_meanarray,HRTFinfo.location_matrix);
   end
   %3. ILDAlone
   
      [Neuron_timelim1{cell_num}.ia_meansurf,...
            Neuron_timelim1{cell_num}.ia_stdsurf,...
            dim2vals,...
            dim1vals,...
            Neuron_timelim1{cell_num}.ia_testpars,...
            spont_spikes, spont_dur, nincl_reps,...
            Neuron_timelim1{cell_num}.ia_locs,...
            Neuron_timelim1{cell_num}.ia_header] = ...
      proc_test899(bird_number,...
      side_of_brain,...
      test_numbers(cell_num,8),...
      1,...
      0,...
      time_limits(time_limit,1),...
      time_limits(time_limit,2) + 20);
   
      [Neuron_timelim1{cell_num}.ia_azi,Neuron_timelim1{cell_num}.ia_ele,temp] = ...
         array2diamond(Neuron_timelim1{cell_num}.ia_meansurf,Neuron_timelim1{cell_num}.ia_locs');
      %Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Neuron_timelim1{cell_num}.ia_azi, Neuron_timelim1{cell_num}.ia_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(Neuron_timelim1{cell_num}.ia_locs(:,1)==el & Neuron_timelim1{cell_num}.ia_locs(:,2)==az))
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
	
   Neuron_timelim1{cell_num}.ia_diamond = intval + temp;
   
   
   %Plot
   
   switch time_limit
      
   case 1,
      
   %Tonal ILD/freq - full time span
   subplot('Position',[0.07 0.85 0.4 0.09])
   set(gca,'FontSize',8);
   plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
   hold on
   colormap(cmap); colorbar
   xlim([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)]);
   ylim([min(Neuron{cell_num}.tif_ildaxis) max(Neuron{cell_num}.tif_ildaxis)]);
   plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'b','LineWidth',2);
   xlabel('Freq (Hz)'); ylabel('ILD (dB)');
   title(['Measured Tonal ILD/freq for cell # ' num2str(cell_num)]);
   
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));

   %BP ILD/freq - full time span
   subplot('Position',[0.55 0.85 0.4 0.09])
   set(gca,'FontSize',8);
   plotsurf(Neuron{cell_num}.bif_freqaxis,Neuron{cell_num}.bif_ildaxis,Neuron{cell_num}.bif_meansurf');
   hold on
   colormap(cmap); colorbar
   xlim([min(Neuron{cell_num}.bif_freqaxis) max(Neuron{cell_num}.bif_freqaxis)]);
   ylim([min(Neuron{cell_num}.bif_ildaxis) max(Neuron{cell_num}.bif_ildaxis)]);
   plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'b','LineWidth',2);
   xlabel('Freq (Hz)'); ylabel('ILD (dB)');
   title(['Measured BP ILD/freq for cell # ' num2str(cell_num)]);
   end
   
   end %end switch
   
   %Tonal ILD/freq - restricted time span
   subplot('Position',[0.07 (1-.15-(time_limit*.15)) 0.4 0.09])
   set(gca,'FontSize',8);
   plotsurf(Neuron_timelim1{cell_num}.tif_freqaxis,...
      Neuron_timelim1{cell_num}.tif_ildaxis,...
      Neuron_timelim1{cell_num}.tif_meansurf');
   hold on
   colormap(cmap); colorbar
   xlim([min(Neuron_timelim1{cell_num}.tif_freqaxis) max(Neuron_timelim1{cell_num}.tif_freqaxis)]);
   ylim([min(Neuron_timelim1{cell_num}.tif_ildaxis) max(Neuron_timelim1{cell_num}.tif_ildaxis)]);
   plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'b','LineWidth',2);
   xlabel('Freq (Hz)'); ylabel('ILD (dB)');
   title(['Measured Tonal ILD/freq, time limits = (' ...
         num2str(Neuron_timelim1{cell_num}.time_limits(1)) ...
         ',' ...
         num2str(Neuron_timelim1{cell_num}.time_limits(2)) ...
         ') ms']);
   
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
   %BP ILD/freq - restricted time span
   subplot('Position',[0.55 (1-.15-(time_limit*.15)) 0.4 0.09])
   set(gca,'FontSize',8);
   plotsurf(Neuron_timelim1{cell_num}.bif_freqaxis,...
      Neuron_timelim1{cell_num}.bif_ildaxis,...
      Neuron_timelim1{cell_num}.bif_meansurf');
   hold on
   colormap(cmap); colorbar
   xlim([min(Neuron_timelim1{cell_num}.bif_freqaxis) max(Neuron_timelim1{cell_num}.bif_freqaxis)]);
   ylim([min(Neuron_timelim1{cell_num}.bif_ildaxis) max(Neuron_timelim1{cell_num}.bif_ildaxis)]);
   plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'b','LineWidth',2);
   xlabel('Freq (Hz)'); ylabel('ILD (dB)');
   title(['Measured BP ILD/freq, time limits = (' ...
         num2str(Neuron_timelim1{cell_num}.time_limits(1)) ...
         ',' ...
         num2str(Neuron_timelim1{cell_num}.time_limits(2)) ...
         ') ms']);   
   end
   end %end loop over time limits
   
   end
   
end %end loop
   
   
   

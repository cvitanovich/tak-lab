clear;close all

bird_number = 899;
side_of_brain = 'r';
%test_numbers = [76 77 79]; %[Tonal_IF BP_IF ILDAlone]
test_numbers = [
   66 68 70 71 65 -57
   76 77 79 80 75 -62
   92 93 94 95 91 -72
   98 100 101 -1 97 -80
   105 106 107 -1 104 -72
   110 111 112 113 109 -75
   123 124 125 126 122 -75
   129 130 131 132 127 -65
   135 136 137 138 134 -73
   143 144 145 146 142 -73
   162 158 159 161 153 -75
   167 168 169 170 166 -70
   173 174 175 176 172 -75]; 

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file = 'out9be';
get_hrtf = 0;

%Get information from bird's HRTF catalogue - for use in predicting ILD Alone Surfaces
   %Header information
   [filetype,...
         info_blocks,...
         n_channels,n_lines, sampling_rate,...
         first_line,last_line,num_locations,...
         comment1, comment2] = mtlrh([hrtf_directory hrtf_file]);
   %Location information
   temp = mtlrdir([hrtf_directory hrtf_file]);
   location_matrix = sph2dbl(temp);
   clear temp
   
   hrtf_freqs = [first_line:(last_line-first_line)/(n_lines-1):last_line]';
   
if(get_hrtf == 1)
   %Transfer Functions
   for loc = 1:num_locations
      temp = mtlrch([hrtf_directory hrtf_file],2*loc-1);
      left_raw(:,loc) = temp;
      left = 20*log10(abs(temp));
      clear temp
      temp = mtlrch([hrtf_directory hrtf_file],2*loc);
      right_raw(:,loc) = temp;
      right = 20*log10(abs(temp));
      clear temp	
      ILD_matrix(:,loc) = right - left;
      ABI_matrix(:,loc) = right + left ./2;
         
      if(mod(loc,10) == 0)
         disp(['Finished location ' num2str(loc)])
      end
   end
   eval(['save ' hrtf_directory num2str(bird_number) 'ildmatrix ILD_matrix left_raw right_raw']);
   eval(['save ' hrtf_directory num2str(bird_number) 'abimatrix ABI_matrix']);
else
   eval(['load ' hrtf_directory num2str(bird_number) 'ildmatrix.mat']);
   eval(['load ' hrtf_directory num2str(bird_number) 'abimatrix.mat']);
end


for cell = 1:size(test_numbers,1);
[Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell,1), 1, 0);


%BP ILD/Freq surface
%Get the mean response surface of the BP ILD/Freq data
[BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell,2), 1, 0);

%Optimal ILD/Freq surface
[pred_ILDAlone,...
      pred_ILDAlone_diam,...
      pred_ILDAlone_azi,...
      pred_ILDAlone_ele,...
      optimal_ildf_surface,faxis,ildaxis] = optimal_ildf(bird_number,side_of_brain,...
   [test_numbers(cell,1) test_numbers(cell,3)],hrtf_file,get_hrtf);

%Get the mean response surface of the ILD Alone data
[ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell,3), 1, 0);
[ILDAlone_azi,ILDAlone_ele,temp] = array2diamond(ILDAlone_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(ILDAlone_azi, ILDAlone_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(locs(:,1)==el & locs(:,2)==az))
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
	
   ILDAlone_diamond = intval + temp;
   ILDAlone_diamond = ILDAlone_diamond./max(max(ILDAlone_diamond));
   
   
   
%Generate the PREDICTED ILD Alone surface from the Optimal Data
[Opt_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
   hrtf_freqs, location_matrix, optimal_ildf_surface, Tonal_dim1vals, Tonal_dim2vals);
[Opt_azi,Opt_ele,Opt_diamond] = array2diamond(Opt_prediction,location_matrix);
Opt_diamond = Opt_diamond./max(max(Opt_diamond));

figure

%Tonal ILD/Freq surface
tonalsurf_h = subplot(3,1,1);
%set(tonalsurf_h,'Position',[0.05 0.05 0.8 0.3]);
plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf');
colorbar
title([num2str(bird_number) side_of_brain ' Depth= ' num2str(Tonal_testpars(7)) ' ITD= '...
      num2str(Tonal_testpars(2)) ' ABI= ' num2str(Tonal_testpars(4))])

BPsurf_h = subplot(3,1,2);
%set(BPsurf_h,'Position',[0.55 0.86 0.3 0.1]);
plotsurf(BP_dim1vals,BP_dim2vals,BP_meansurf');
colorbar

optsurf_h = subplot(3,1,3);
%set(optsurf_h,'Position',[0.05 0.55 0.8 0.3]);
plotsurf(faxis,ildaxis,optimal_ildf_surface');
colorbar
xlabel('Frequency (Hz)');
ylabel('ILD (dB)');


figure
subplot(1,3,1)
plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diamond);
title('Measured')
subplot(1,3,2)
plotdiam(Opt_azi,Opt_ele,Opt_diamond);
title('Predicted (optimal)')
subplot(1,3,3)
plotdiam(pred_ILDAlone_azi,pred_ILDAlone_ele,pred_ILDAlone_diam);
title('Predicted (reverse)')



end %end for loop


%Script to test various ways of producing an "optimal" ILD/Freq response surface from
%a measured ILDAlone spatial response surface

clear;close all
nntwarn off %needed to suppress warning messages regarding obsolete neural network functions

%Input needed for functions below
bird_number = 899;
side_of_brain = 'r';
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


%Begin loop over cells
for cell = 1:size(test_numbers,1)
%for cell = 1:5
   %1. Get the cell's Tonal ILD/Freq RS
   [Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,1), 1, 0);
   
   %2. Get the cell's BP ILD/Freq RS
   [BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,2), 1, 0);
   
   %3. Get the cell's "optimal" ILD/Freq RS using the method of Weighted Mean
   %this also produces the "reverse" of the pure line integral (pred_ILDAlone series)
   [pred_ILDAlone,...
         pred_ILDAlone_diam,...
         pred_ILDAlone_azi,...
         pred_ILDAlone_ele,...
         optimal_ildf_surface,faxis,ildaxis] = optimal_ildf(bird_number,side_of_brain,...
      [test_numbers(cell,1) test_numbers(cell,3)],hrtf_file,get_hrtf);
   

   %4. Get the mean response surface of the ILD Alone data
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
   
   %5. Generate the PREDICTED ILDAlone RS from the Weighted Mean "optimal" surface
   [Opt_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, optimal_ildf_surface, Tonal_dim1vals, Tonal_dim2vals);
   [Opt_azi,Opt_ele,Opt_diamond] = array2diamond(Opt_prediction,location_matrix);
   Opt_diamond = Opt_diamond./max(max(Opt_diamond));
   
   %6. Generate the PREDICTED ILDAlone RS from the BackPropagation "optimal" surface
   figure(1) %needed for the plotting that takes place by the nnet training
   [bpopt_ILDAlone_azi,...
      bpopt_ILDAlone_ele,...
      bpopt_ILDAlone_diamond,...
      bp_optimal_ildf_surface] = ...
   bp_optimal_ildf(bird_number,side_of_brain,[test_numbers(cell,1) test_numbers(cell,3)],...
   hrtf_file,get_hrtf);

   [bpopt_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, bp_optimal_ildf_surface, Tonal_dim1vals, Tonal_dim2vals);
   [bpopt_azi,bpopt_ele,bpopt_diamond] = array2diamond(bpopt_prediction,location_matrix);
   bpopt_diamond = bpopt_diamond./max(max(bpopt_diamond));


   %Plotting:
   figure %ILD/Freq surfaces
   
   %Tonal ILD/Freq surface
   tonalsurf_h = subplot(4,1,1);
   %set(tonalsurf_h,'Position',[0.05 0.05 0.8 0.3]);
   plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf');
   colorbar
   title([num2str(bird_number) side_of_brain ' Depth= ' num2str(Tonal_testpars(7)) ' ITD= '...
         num2str(Tonal_testpars(2)) ' ABI= ' num2str(Tonal_testpars(4))])
   
   %BP ILD/Freq surface
   BPsurf_h = subplot(4,1,2);
   %set(BPsurf_h,'Position',[0.55 0.86 0.3 0.1]);
   plotsurf(BP_dim1vals,BP_dim2vals,BP_meansurf');
   colorbar
   
   %Weighted Mean Optimal ILD/Freq surface
   optsurf_h = subplot(4,1,3);
   %set(optsurf_h,'Position',[0.05 0.55 0.8 0.3]);
   plotsurf(faxis,ildaxis,optimal_ildf_surface');
   colorbar
   
   %BackProp Optimal ILD/Freq surface
   bp_optsurf_h = subplot(4,1,4);
   plotsurf(Tonal_dim1vals,Tonal_dim2vals,bp_optimal_ildf_surface');
   colorbar
   xlabel('Frequency (Hz)');
   ylabel('ILD (dB)');
   
   figure %ILDAlone RS's
   subplot(2,2,1)
   plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diamond);
   title('Measured')
   
   for cc_loc = 1:length(locs)
      cc_index(cc_loc) = max(find(location_matrix(1,:) == locs(cc_loc,1)...
         & location_matrix(2,:) == locs(cc_loc,2)));
   end
   
   subplot(2,2,2)
   plotdiam(pred_ILDAlone_azi,pred_ILDAlone_ele,pred_ILDAlone_diam);
   cc = corrcoef(ILDAlone_meansurf,pred_ILDAlone);
   cc = cc(1,2);
   title(['Pure Line Integrals, r = ' num2str(cc)])
   
   subplot(2,2,3)
   plotdiam(Opt_azi,Opt_ele,Opt_diamond);
   cc = corrcoef(ILDAlone_meansurf,Opt_prediction(cc_index));
   cc = cc(1,2);
   title(['Weighted Mean Prediction, r = ' num2str(cc)])
   
   subplot(2,2,4)
   plotdiam(bpopt_azi,bpopt_ele,bpopt_diamond);
   cc = corrcoef(ILDAlone_meansurf,bpopt_prediction(cc_index));
   cc = cc(1,2);
   title(['Backprop Prediction, r = ' num2str(cc)])
   

end %end for loop over cells


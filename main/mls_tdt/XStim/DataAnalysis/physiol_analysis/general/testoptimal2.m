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
colormap_var = 'hot';
plotflag = 0;

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

faxis = 2000:475:11500;
ildaxis = -30:4:30;
[XI,YI] = meshgrid(faxis,ildaxis);
%Begin loop over cells
for cell = 1:size(test_numbers,1)
%for cell = 1:5
   %1. Get the cell's Tonal ILD/Freq RS
   [Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,1), 1, 0);
   Tonal_meansurf = Tonal_meansurf./max(max(abs(Tonal_meansurf)));
   if(size(Tonal_meansurf,1) ~= length(faxis) | size(Tonal_meansurf,2) ~= length(ildaxis))
      Tonal_ildf(:,:,cell) = interp2(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf',XI',YI');
   else
      Tonal_ildf(:,:,cell) = Tonal_meansurf;
   end
   %Tonal_meansurf = Tonal_meansurf - mean(dezero(Tonal_meansurf(:,1))); %estimate and subtract spontaneous rate
   
   %2. Generate the PREDICTED ILD Alone surface from the Tonal Data
   [Tonal_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, Tonal_meansurf, Tonal_dim1vals, Tonal_dim2vals);
   [Tonal_azi,Tonal_ele,Tonal_diamond] = array2diamond(Tonal_prediction,location_matrix);
   Tonal_diamond = Tonal_diamond./max(max(Tonal_diamond));
   
   %3. Get the cell's BP ILD/Freq RS
   [BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,2), 1, 0);
   BP_meansurf = BP_meansurf./max(max(abs(BP_meansurf)));
   if(size(BP_meansurf,1) ~= length(faxis) | size(BP_meansurf,2) ~= length(ildaxis))
      BP_ildf(:,:,cell) = interp2(BP_dim1vals,BP_dim2vals,BP_meansurf',XI',YI');
   else
      BP_ildf(:,:,cell) = BP_meansurf;
   end
   %BP_meansurf = BP_meansurf - mean(dezero(BP_meansurf(:,1))); %estimate and subtract spontaneous rate
   
   %4. Generate the PREDICTED ILD Alone surface from the BP Data
   [BP_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, BP_meansurf, BP_dim1vals, BP_dim2vals);
   [BP_azi,BP_ele,BP_diamond] = array2diamond(BP_prediction,location_matrix);
   BP_diamond = BP_diamond./max(max(BP_diamond));
   
   %5. Get the mean response surface of the ILD Alone data
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
   
   %6. Generate the PREDICTED ILDAlone RS from the BackPropagation "optimal" surface
   figure(1) %needed for the plotting that takes place by the nnet training
   [bpopt_ILDAlone_azi,...
      bpopt_ILDAlone_ele,...
      bpopt_ILDAlone_diamond,...
      bp_optimal_ildf_surface] = ...
   bp_optimal_ildf(bird_number,side_of_brain,[test_numbers(cell,1) test_numbers(cell,3)],...
   hrtf_file,get_hrtf);
   bp_optimal_ildf_surface = bp_optimal_ildf_surface./max(max(abs(bp_optimal_ildf_surface)));
   if(size(bp_optimal_ildf_surface,1) ~= length(faxis) | size(bp_optimal_ildf_surface,2) ~= length(ildaxis))
      BackProp_ildf(:,:,cell) = interp2(Tonal_dim1vals,Tonal_dim2vals,bp_optimal_ildf_surface',XI',YI');
   else
      BackProp_ildf(:,:,cell) = bp_optimal_ildf_surface;
   end

   [bpopt_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, bp_optimal_ildf_surface, Tonal_dim1vals, Tonal_dim2vals);
   [bpopt_azi,bpopt_ele,bpopt_diamond] = array2diamond(bpopt_prediction,location_matrix);
   bpopt_diamond = bpopt_diamond./max(max(bpopt_diamond));

   if (plotflag == 1)
   %Plotting:
   figure %ILD/Freq surfaces
   
   %Tonal ILD/Freq surface
   tonalsurf_h = subplot(3,1,1);
   %set(tonalsurf_h,'Position',[0.05 0.05 0.8 0.3]);
   plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf');
   colormap(colormap_var);
   colorbar
   title([num2str(bird_number) side_of_brain ' Depth= ' num2str(Tonal_testpars(7)) ' ITD= '...
         num2str(Tonal_testpars(2)) ' ABI= ' num2str(Tonal_testpars(4))])
   text(2000,20,'\fontname{times}Tonal ILD/Freq RS','FontSize',10,'Color','white')
   
   %BP ILD/Freq surface
   BPsurf_h = subplot(3,1,2);
   %set(BPsurf_h,'Position',[0.55 0.86 0.3 0.1]);
   plotsurf(BP_dim1vals,BP_dim2vals,BP_meansurf');
   colormap(colormap_var);
   colorbar
   text(2000,20,'\fontname{times}BandPass ILD/Freq RS','FontSize',10,'Color','white')
   
   %BackProp Optimal ILD/Freq surface
   bp_optsurf_h = subplot(3,1,3);
   plotsurf(Tonal_dim1vals,Tonal_dim2vals,bp_optimal_ildf_surface');
   colormap(colormap_var);
   colorbar
   text(2000,20,'\fontname{times}BackProp ILD/Freq RS','FontSize',10,'Color','white')
   xlabel('Frequency (Hz)');
   ylabel('ILD (dB)');
   
   figure %ILDAlone RS's
   
   %Measured
   subplot(2,2,1)
   plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diamond);
   colormap(colormap_var);
   title('Measured')
   
   for cc_loc = 1:length(locs)
      cc_index(cc_loc) = max(find(location_matrix(1,:) == locs(cc_loc,1)...
         & location_matrix(2,:) == locs(cc_loc,2)));
   end
   
   %Tonal
   subplot(2,2,2)
   plotdiam(Tonal_azi,Tonal_ele,Tonal_diamond);
   colormap(colormap_var);
   cc = corrcoef(ILDAlone_meansurf,Tonal_prediction(cc_index));
   cc = cc(1,2);
   title(['Tonal Prediction, r = ' num2str(cc)])
   
   %Bandpassed
   subplot(2,2,3)
   plotdiam(BP_azi,BP_ele,BP_diamond);
   colormap(colormap_var);
   cc = corrcoef(ILDAlone_meansurf,BP_prediction(cc_index));
   cc = cc(1,2);
   title(['BandPassed Prediction, r = ' num2str(cc)])
   
   %BackProp
   subplot(2,2,4)
   plotdiam(bpopt_azi,bpopt_ele,bpopt_diamond);
   colormap(colormap_var);
   cc = corrcoef(ILDAlone_meansurf,bpopt_prediction(cc_index));
   cc = cc(1,2);
   title(['Backprop Prediction, r = ' num2str(cc)])
   
   end %plotflag
   

end %end for loop over cells

ildf_diff_tonal = BackProp_ildf(1:size(Tonal_ildf,1),1:size(Tonal_ildf,2),:) - Tonal_ildf;
ildf_diff_tonal(find(isnan(ildf_diff_tonal))) = 0;
mean_ildf_diff_tonal = mean(ildf_diff_tonal,3);
figure
subplot(2,1,1)
plotsurf(Tonal_dim1vals,Tonal_dim2vals,mean_ildf_diff_tonal');
colormap(colormap_var);
colorbar
text(2000,20,'\fontname{times}Tonal ILD/Freq Diff','FontSize',10,'Color','white')



ildf_diff_bp = BackProp_ildf(1:size(BP_ildf,1),1:size(BP_ildf,2),:) - BP_ildf;
ildf_diff_bp(find(isnan(ildf_diff_bp))) = 0;
mean_ildf_diff_bp = mean(ildf_diff_bp,3);
subplot(2,1,2)
plotsurf(BP_dim1vals,BP_dim2vals,mean_ildf_diff_bp');
colormap(colormap_var);
colorbar
text(2000,20,'\fontname{times}BP ILD/Freq Diff','FontSize',10,'Color','white')


return

for fig = 4:27
   figure(fig)
   print
end


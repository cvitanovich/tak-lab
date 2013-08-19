%Script to study whether the ILD-distance method is *predictive*
%That is, can a subset of locations be used to generate Activity vs. d_ILD curves that predict
%the activity for the rest of the locations?

%Load neuronal data
if(~exist('Neuron','var'))
   load 'd:\mlspezio\matlab\save\Neuron_39_regress_twelfth'
end
Neuron = Neuron_twelfth;

%Load HRTF files for 901, 913 & 894
if(exist('d:\mlspezio\matlab\save\HRTFAll_3.mat'))
   eval(['load d:\mlspezio\matlab\save\HRTFAll_3'])
end

vert_pos = [0.85 0.68 0.31 0.03];

frac_points = 0.3; %Choose what fraction of the # of total points to include
%for parametrizing the model
for cell_num = [22:39]
   
   clear ILD_matrix HRTFinfo ILDf_freqaxis ILDf_ildaxis ILDf_meansurf ILDf_stdsurf ILD_matrix_focus
   clear IA_meansurf IA_stdsurf IA_diamond IA_locs IA_azi IA_ele bbif_meansurf
   
   %Specify Site Parameters
   if(cell_num < 25)
      bird_number = 899;
      depth = Neuron{cell_num}.tif_testpars(7);
      site_number = cell_num;
      side_of_brain = '';
      stim_dur = Neuron{cell_num}.tif_testpars(8);
      isi = [];
      session_num = [];
   else
      bird_number = Neuron{cell_num}.abi_testpars.bird_number;
      session_num = Neuron{cell_num}.abi_testpars.session_num;
      depth = Neuron{cell_num}.abi_testpars.depth;
      site_number = Neuron{cell_num}.abi_testpars.recording_site;
      side_of_brain = Neuron{cell_num}.abi_testpars.side_of_brain;
      stim_dur = Neuron{cell_num}.abi_testpars.curr_stimdur;
      isi = Neuron{cell_num}.abi_testpars.test_ISI;
   end
   
   %Choose HRTF file
   if(cell_num < 25) %for 899
      if(0)
         hrtf_file = 'd:\mlspezio\owl_hrtfdata\899\out9be';
         get_HRTF = 1;
         get_ITDmatrix = 0;
         hrtf_file2 = [];
         [ILD_matrix,ABL_matrix,ITD_matrix,HRTFinfo] = ...
            get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
         save 'd:\mlspezio\matlab\save\899HRTF' ILD_matrix ABL_matrix HRTFinfo
      else
         load 'd:\mlspezio\matlab\save\899HRTF' ILD_matrix ABL_matrix HRTFinfo
      end
   else
      switch bird_number
      case 901
         ILD_matrix = HRTFAll{1}.ILD_matrix;
         HRTFinfo = HRTFAll{1}.HRTFinfo;
      case 913
         ILD_matrix = HRTFAll{2}.ILD_matrix;
         HRTFinfo = HRTFAll{2}.HRTFinfo;
      case 894
         ILD_matrix = HRTFAll{3}.ILD_matrix;
         HRTFinfo = HRTFAll{3}.HRTFinfo;
      end %end switch
   end %end if
   
   %Calculate ILD-distance-based measures
   disp(['Processing cell # ' num2str(cell_num)])
   
   %Choose the elements for analysis:
   if(isfield(Neuron{cell_num},'gif_freqaxis'))
      ILDf_freqaxis  = Neuron{cell_num}.gif_freqaxis;
      ILDf_ildaxis   = Neuron{cell_num}.gif_ildaxis;
      ILDf_meansurf  = Neuron{cell_num}.gif_mean;
      ILDf_meanarray = Neuron{cell_num}.gia_meanarray;
   elseif(isfield(Neuron{cell_num},'tif_freqaxis'))
      ILDf_freqaxis  = Neuron{cell_num}.tif_freqaxis;
      ILDf_ildaxis   = Neuron{cell_num}.tif_ildaxis;
      ILDf_meansurf  = Neuron{cell_num}.tif_meansurf;
      ILDf_meanarray = Neuron{cell_num}.tia_meanarray;
   end
   
   if(cell_num < 25)
      IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
      IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
      IA_diamond  = Neuron{cell_num}.ia_diamond{1};
      IA_locs		= Neuron{cell_num}.ia_locs{1};
      IA_azi		= Neuron{cell_num}.ia_azi{1};
      IA_ele		= Neuron{cell_num}.ia_ele{1};
   else
      IA_meansurf = Neuron{cell_num}.ia_mean{1};
      IA_stdsurf  = Neuron{cell_num}.ia_std{1};
      IA_diamond  = Neuron{cell_num}.ia_diamond{1};
      IA_locs		= Neuron{cell_num}.ia_locs{1};
      IA_azi		= Neuron{cell_num}.ia_azi{1};
      IA_ele		= Neuron{cell_num}.ia_ele{1};
   end
   
   if(size(IA_locs,1) < 3)
      IA_locs = IA_locs';
   end
   
   %Get those locations in the ILDmatrix that match the measured ILDAlone locations
   for num_loc = 1:size(IA_locs,1)
      ILDmat_index(num_loc) =...
         max(find(HRTFinfo.location_matrix(1,:) == IA_locs(num_loc,1) &...
         HRTFinfo.location_matrix(2,:) == IA_locs(num_loc,2)));
   end
   ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   
   %Get the minimum and maximum frequency indices
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      ILDf_freqaxis(length(ILDf_freqaxis))));
   
   if(isfield(Neuron{cell_num},'ts_meansurf'))
      TS_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ts_meansurf{1} == max(Neuron{cell_num}.ts_meansurf{1})));
   elseif(isfield(Neuron{cell_num},'ts_mean'))
      TS_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ts_mean{1} == max(Neuron{cell_num}.ts_mean{1})));
   end
   IA_bestloc_ind(cell_num) = max(find(IA_meansurf == max(IA_meansurf)));
   if(cell_num == 1) IA_bestloc_ind(cell_num) = 129; end
   if(cell_num == 4) IA_bestloc_ind(cell_num) = 198; end
   if(cell_num == 7) IA_bestloc_ind(cell_num) = 163; end
   if(cell_num == 9) IA_bestloc_ind(cell_num) = 282; end
   if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
   if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
   if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
   if(cell_num == 18) IA_bestloc_ind(cell_num) = 78; end
   if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
   
   %Determine Frequency Bands to use in the MR analysis
   clear freq_bands
   for freq_num = 1:length(ILDf_freqaxis)
      [freq_bands(freq_num,1) freq_bands(freq_num,2)] = bandlimits(ILDf_freqaxis(freq_num),1/12);
      if(freq_bands(freq_num,1) < 2000) freq_bands(freq_num,1) = 2000; end
      if(freq_bands(freq_num,2) < 2000) freq_bands(freq_num,2) = 2000; end
      if(freq_bands(freq_num,1) > 11000) freq_bands(freq_num,1) = 11000; end
      if(freq_bands(freq_num,2) > 11000) freq_bands(freq_num,2) = 11000; end
   end
   
   %Calculate distance_ILD matrix, unsigned and signed
   disp('Calculating signed distance matrix...')
   sign_flag = 1;
   dist_ild_signed = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
      IA_bestloc_ind(cell_num),freq_bands,sign_flag);
   [dist_ild_s_sort,dsort_s_index] = sort(dist_ild_signed,2);
   %***Important: reverse the sorting for the comparison
   [junk,reverse_sort_index] = sort(dsort_s_index,2);
   clear junk
   
   %Get the indices for inclusion in the Activity vs. d_ILD calculation
   randtemp = randperm(size(ILD_matrix_focus,2));
   keep_inds = randtemp(1:round(frac_points*size(ILD_matrix_focus,2)));
   ind_num_in = 0;
   ind_num_out = 0;
   for loc = 1:size(IA_meansurf,2)
      if(~all(loc - keep_inds))
         ind_num_in = ind_num_in + 1;
         train_inds(ind_num_in) = loc;
      else
         ind_num_out = ind_num_out + 1;
         test_inds(ind_num_out) = loc;
      end
   end
   
   %Calculate Activity vs. ILD_distance curves
   alpha = 0.15;
   for freq_num = 1:length(ILDf_freqaxis)
      %Determine how to sort the training indices, and how to reverse the sort
      [y,key_index] = sort(dist_ild_signed(freq_num,train_inds));
      [y,rev_key_index] = sort(key_index);
      %Find max and min d_ILDs, and use them to get the curve
      [y,min_ind] = min(dist_ild_signed(freq_num,:));
      [y,max_ind] = max(dist_ild_signed(freq_num,:));
      if(all(train_inds - min_ind) & all(train_inds - max_ind))
         [actd_curve,resia] = ...
            locregress1([dist_ild_signed(freq_num,min_ind) ...
               sort(dist_ild_signed(freq_num,train_inds)) ...
               dist_ild_signed(freq_num,max_ind)],...
            IA_meansurf([min_ind train_inds(key_index) max_ind]),alpha);
         temp_curve = actd_curve(2:length(actd_curve)-1);
         d_ild_all{freq_num} = [dist_ild_signed(freq_num,min_ind) ...
               sort(dist_ild_signed(freq_num,train_inds)) ...
               dist_ild_signed(freq_num,max_ind)];
         actd_curve_all{freq_num} = actd_curve;
         actd_curve_use(freq_num,:) = temp_curve(rev_key_index);
      elseif(all(train_inds - min_ind))
         [actd_curve,resia] = ...
            locregress1([dist_ild_signed(freq_num,min_ind) ...
               sort(dist_ild_signed(freq_num,train_inds))],...
            IA_meansurf([min_ind train_inds(key_index)]),alpha);
         temp_curve = actd_curve(2:length(actd_curve));
         d_ild_all{freq_num} = [dist_ild_signed(freq_num,min_ind) ...
               sort(dist_ild_signed(freq_num,train_inds))];
         actd_curve_all{freq_num} = actd_curve;
         actd_curve_use(freq_num,:) = temp_curve(rev_key_index);
      elseif(all(train_inds - max_ind))
            [actd_curve,resia] = ...
               locregress1([sort(dist_ild_signed(freq_num,train_inds)) ...
               dist_ild_signed(freq_num,max_ind)],...
               IA_meansurf([train_inds(key_index) max_ind]),alpha);
            temp_curve = actd_curve(1:length(actd_curve)-1);
            d_ild_all{freq_num} = [sort(dist_ild_signed(freq_num,train_inds)) ...
               dist_ild_signed(freq_num,max_ind)];
            actd_curve_all{freq_num} = actd_curve;
            actd_curve_use(freq_num,:) = temp_curve(rev_key_index);
      else
            [actd_curve,resia] = ...
               locregress1(sort(dist_ild_signed(freq_num,train_inds)),...
               IA_meansurf(train_inds(key_index)),alpha);
            d_ild_all{freq_num} = sort(dist_ild_signed(freq_num,train_inds));
            actd_curve_all{freq_num} = actd_curve;
            actd_curve_use(freq_num,:) = actd_curve(rev_key_index);
      end
         disp(['Finished bb simulation for ' num2str(ILDf_freqaxis(freq_num)) ' Hz'])
   end%end loop over frequencies
      
      %Perform multiple regression - stepwise or direct
      alpha = 0.05; %confidence interval criterion
      
      if(0) %run a direct regression
         [beta,betaci,r,rint,tvals,stats] = regress(IA_meansurf(train_inds_sort)',...
            actd_curve_use',alpha);
      else %run a stepwise regression
         stepwise(actd_curve_use',...
            IA_meansurf(train_inds)',...
            1:size(freq_bands,1),alpha);
         keyboard
      end
      
      %Only take those coefficients that contribute significantly
      newbeta = beta;
      for freq_num = 1:length(ILDf_freqaxis)
         if( (sign(betaci(freq_num,1)) ~= sign(betaci(freq_num,2))) | ~all(out - freq_num)) %out is a vector of freq_bands not included in the model
            newbeta(freq_num) = 0;
         end
      end
      
   %Test agreement on remaining locations
   for freq_num = 1:length(ILDf_freqaxis)
      activities(freq_num,:) = interp1(d_ild_all{freq_num},actd_curve_all{freq_num},dist_ild_signed(freq_num,:));
   end
   
   IA_predicted = newbeta' * activities;
   
   temp = corrcoef(IA_meansurf(test_inds),IA_predicted(test_inds));
   Perc_30_var(cell_num) = temp(1,2)^2;
   
   %Compute multicollinearity of actd_curve_use
   for m=1:size(actd_curve_use,1)
      for n=1:size(actd_curve_use,1)
         temp = corrcoef(actd_curve_use(m,:),actd_curve_use(n,:));
         cc{cell_num}(m,n) = temp(1,2);
      end
   end
   
end %end loop over cells

for cell_num = [1:2 4:11 13:20 22:24]
   %Calculate ILD-distance-based measures
   disp(['Processing cell # ' num2str(cell_num)])
   
   %Choose the elements for analysis:
   if(isfield(Neuron{cell_num},'gif_freqaxis'))
      ILDf_freqaxis  = Neuron{cell_num}.gif_freqaxis;
   elseif(isfield(Neuron{cell_num},'tif_freqaxis'))
      ILDf_freqaxis  = Neuron{cell_num}.tif_freqaxis;
   end
   
   figure
   pcolor(ILDf_freqaxis,ILDf_freqaxis,cc{cell_num});
   colorbar
   title(['Cell Num = ' num2str(cell_num)])
   colormap(1-gray)
   print
   pause(2)
end
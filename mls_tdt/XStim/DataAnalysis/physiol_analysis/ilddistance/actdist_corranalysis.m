%Linmodel9 - script descibing and testing a linear model of spectral integration in ILDAlone RS
clear;
Model = 'Linmodel9_{LR}';

%set parameters
bird_number = 899;
begin_cells = 1;
end_cells = 24;

%Graphics parameters
cmap = jet;
linecolor = 'white';
bkgndcolor = [0.85 0.85 0.85];
fig_ht = 10;
fig_wd = 8;

%Printing
print_flag = 0;

%Directories & files
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\HRTF_files\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;

%Local Regression Parameters
alpha_dist = 0.05; alpha_ILD = 0.05; mode = 0;

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix test_numbers reps

for cell_num = begin_cells:end_cells
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
      
      %Choose the elements for analysis:
      %1. ILD_frequency frequency axis, 2. ILDAlone meansurf and diamond
      
      clear ILDf_freqaxis cc_actfreq
      if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
         ILDf_freqaxis = Neuron{cell_num}.bif_freqaxis;
      else
         ILDf_freqaxis = Neuron{cell_num}.tif_freqaxis;
      end
      
         IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
         IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
         IA_diamond  = Neuron{cell_num}.ia_diamond{1};
         IA_locs		= Neuron{cell_num}.ia_locs{1};
         IA_azi		= Neuron{cell_num}.ia_azi{1};
         IA_ele		= Neuron{cell_num}.ia_ele{1};
      
   disp(['Processing cell # ' num2str(cell_num)])
   %Get the contour for the SRF of the neuron
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
   [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
      srfstats(Neuron{cell_num}.ts_azi{1},...
      Neuron{cell_num}.ts_ele{1},...
      Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
   end
   
   %Get those locations in the ILDmatrix that match the measured ILDAlone locations
   for num_loc = 1:size(IA_locs,1)
      ILDmat_index(num_loc) =...
         max(find(HRTFinfo.location_matrix(1,:) == IA_locs(num_loc,1) &...
         HRTFinfo.location_matrix(2,:) == IA_locs(num_loc,2)));
   end
   ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   
   %Get the frequency indices
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      ILDf_freqaxis(length(ILDf_freqaxis))));
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
   
   disp(['Modelling activity for ' num2str(length(IA_meansurf)) ' locations...'])
   clear dist_ild dist_ild_signed dist_ild_signed_sort d_s_sort_index result_pieces
   clear residual result ActILD_result ActILD_residual ILD_matrix_focus_restr_sort ILDf_linpredict
   clear IA_linpredict ILD_bands freq_bands
   bandwidth = 1/3; %units of octaves
   
   %Determine Frequency Bands (from Tonal ILD/frequency
   for freq_num = 1:length(ILDf_freqaxis)
      %Determine frequency limits
      [freq_bands(freq_num,1),freq_bands(freq_num,2)] = ...
         bandlimits(ILDf_freqaxis(freq_num),bandwidth);
   end
   
   %Calculate the mean ILD across the frequency bands
   disp('Getting ILD across frequency bands...')
   sign_flag = 1;
   [ILD_bands] = get_ild_bands(ILD_matrix_focus,HRTFinfo.hrtf_freqs,freq_bands,sign_flag);

   %Calculate distance_ILD matrix, unsigned and signed
   disp('Calculating distance matrices...')
   sign_flag = 0;
   dist_ild = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
      IA_bestloc_ind(cell_num),freq_bands,sign_flag);
   sign_flag = 1;
   dist_ild_signed = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
      IA_bestloc_ind(cell_num),freq_bands,sign_flag);
   
   
   for freq_num = 1:length(ILDf_freqaxis)
      [y,freqind] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(freq_num)));
      %Calculate correlation between Activity & ILD_distance
      temp = corrcoef(dist_ild(freq_num,:),IA_meansurf);
      cc_actfreq(freq_num) = -1*temp(1,2);
   end
   
      figure
      maxind = find(cc_actfreq == max(cc_actfreq));
      a = corrcoef(dist_ild');
      norm_cc_actfreq = cc_actfreq/max(cc_actfreq);
      plot(ILDf_freqaxis,norm_cc_actfreq);
      hold on
      plot(ILDf_freqaxis,a(maxind,:),'r--');
      temp1 = norm_cc_actfreq - a(maxind,:);
      temp1(maxind) = 1;
      plot(ILDf_freqaxis,temp1,'g-.')
      
      tonal_freqcurve = max(Neuron{cell_num}.tif_meansurf,[],2);
      tonal_freqcurve = tonal_freqcurve/max(tonal_freqcurve);
      plot(Neuron{cell_num}.tif_freqaxis,tonal_freqcurve,'m')
      title(['Correlations for cell # ' num2str(cell_num)]);
   
end

end

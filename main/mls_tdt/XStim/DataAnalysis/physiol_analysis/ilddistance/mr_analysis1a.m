%mr_analysis1
%Script to carry out a Multiple Regression modelling and analysis of ILDAlone activity
%Dependent variable: ILDAlone activity
%Independent variables:	activity as determined by the Local Regression of Activity on signed ILD_distance

clear; close all;

%set parameters
bird_number = 899;
begin_cells = 24;
tot_num_cells = 1;
end_cells = begin_cells + tot_num_cells - 1;

%Graphics parameters
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

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix test_numbers reps

for cell_num = begin_cells:end_cells
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
      
      disp(['Processing cell # ' num2str(cell_num)])
      
      %Choose the elements for analysis:     
      ILDf_freqaxis  = Neuron{cell_num}.tif_freqaxis;
      ILDf_ildaxis   = Neuron{cell_num}.tif_ildaxis;
      ILDf_meansurf  = Neuron{cell_num}.tif_meansurf;
      ILDf_meanarray = Neuron{cell_num}.tia_meanarray;
      
      IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
      IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
      IA_diamond  = Neuron{cell_num}.ia_diamond{1};
      IA_locs		= Neuron{cell_num}.ia_locs{1};
      IA_azi		= Neuron{cell_num}.ia_azi{1};
      IA_ele		= Neuron{cell_num}.ia_ele{1};
      
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
   for freq_num = 1:length(ILDf_freqaxis)-1
      freq_bands(freq_num,1) = ILDf_freqaxis(freq_num);
      freq_bands(freq_num,2) = ILDf_freqaxis(freq_num + 1);
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
   
   %Calculate Activity vs. ILD_distance curves
   alpha = 0.15;
   for freq_num = 1:length(ILDf_freqaxis)-1
      [tempia,resia] = locregress1(dist_ild_s_sort(freq_num,:),...
         IA_meansurf(dsort_s_index(freq_num,:)),alpha);
      bbia_meansurf(freq_num,:) = tempia(reverse_sort_index(freq_num,:));
      disp(['Finished bb simulation for ' num2str(ILDf_freqaxis(freq_num)) ' Hz'])
   end
   
   %Perform multiple regression
   alpha = 0.01; %confidence interval criterion
   stepwise(bbia_meansurf',IA_meansurf',1:length(freq_bands),alpha);
   [b,bint,r,rint,tvals,stats] = regress(IA_meansurf',bbia_meansurf',alpha);
   
   cctemp = corrcoef(IA_meansurf,ILDf_meanarray(ILDmat_index));
   perc_var = cctemp(1,2)^2;
   
   h = figure;
   set(h,'Units','normal','Position',[0.01,0.01,0.4,0.4]);
   plot(ILDf_freqaxis,max(ILDf_meansurf,[],2),'*-');
   xlabel('Frequency (Hz)'); ylabel('Spikes/stim');
   title(['Var explained by linear combination = ' num2str(perc_var) '%'])
   
   end %end if

end %end for loop

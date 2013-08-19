%Script to test the idea that simple averaging in frequency is enough to explain why
%complex noises are better than simple tones at predicting a cell's ILDAlone RS

clear;close all

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

for cell = 1:size(test_numbers,1)
   %1. Get the cell's Tonal ILD/Freq RS
   [Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,1), 1, 0);
   Tonal_meansurf = Tonal_meansurf./max(max(abs(Tonal_meansurf)));
   Tonal_freqcurve = max(Tonal_meansurf,[],2);
   
   %2. Generate the PREDICTED ILDAlone surface from Tonal Data
   [Tonal_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, Tonal_meansurf, Tonal_dim1vals, Tonal_dim2vals);
   
   %3. Get the cell's BP ILD/Freq RS
   [BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,2), 1, 0);
   BP_meansurf = BP_meansurf./max(max(abs(BP_meansurf)));

   %4. Generate the PREDICTED ILD Alone surface from the BP Data
   [BP_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
      hrtf_freqs, location_matrix, BP_meansurf, BP_dim1vals, BP_dim2vals);
   
   %5. Get the mean response surface of the ILD Alone data
   [ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
      proc_test899(bird_number, side_of_brain, test_numbers(cell,3), 1, 0);
   
   for cc_loc = 1:length(locs)
      cc_index(cc_loc) = max(find(location_matrix(1,:) == locs(cc_loc,1)...
         & location_matrix(2,:) == locs(cc_loc,2)));
   end
   
   count = 1;
   temp = corrcoef(ILDAlone_meansurf,Tonal_prediction(cc_index));
   cc_vec(cell,count) = temp(1,2);
   temp = corrcoef(ILDAlone_meansurf,BP_prediction(cc_index));
   cc_bp_vec(cell) = temp(1,2);
   
   figure
   subplot(11,1,count)
   plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf');
   colormap(colormap_var);
   colorbar

   
   for bpwidth = 0.1:0.1:1 %from 0.1 to 1 octave
      count = count+1;
      %Do the appropriate averaging
      for num_freq = 1:length(Tonal_dim1vals)
         [lo_limit,hi_limit] = bandlimits(Tonal_dim1vals(num_freq),bpwidth);
         index_1 = min(find(Tonal_dim1vals >= lo_limit));
         index_2 = max(find(Tonal_dim1vals <= hi_limit));
         %mult_mat = repmat(Tonal_freqcurve(index_1:index_2),1,size(Tonal_meansurf,2));
         Tonal_meansurf_new(num_freq,:) = mean(Tonal_meansurf(index_1:index_2,:),1);
      end
      
      subplot(11,1,count)
      plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf_new');
      colormap(colormap_var);
      colorbar
      
      %Generate the PREDICTED ILDAlone surface from Tonal Data
      [Tonal_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
         hrtf_freqs, location_matrix, Tonal_meansurf_new, Tonal_dim1vals, Tonal_dim2vals);
      clear Tonal_meansurf_new
      
      temp = corrcoef(ILDAlone_meansurf,Tonal_prediction(cc_index));
      cc_vec(cell,count) = temp(1,2);
      
   end
   
   figure
   hold on
   plot(0:0.1:1,cc_vec(cell,:))
   plot(0.3,cc_bp_vec(cell))
   text(0.3,cc_bp_vec(cell),'\bullet\leftarrow\fontname{times}For measured 0.3 BP',...
      'FontSize',10,'Color','black')
   text(0.3,cc_vec(cell,4),'\leftarrow\fontname{times}Averaged Tonal ILD/Freq',...
      'FontSize',10,'Color','black')
   xlabel('Bandwidth of Averaging (octaves)')
   ylabel('Correlation between ILDAlone & ILD/Freq Prediction');
   
end

   figure
   hold on
   for cell = 1:size(test_numbers,1)
      plot(0:0.1:1,cc_vec(cell,:),'b.')
      plot(0.3,cc_bp_vec(cell),'r*')
   end
   xlabel('Bandwidth of Averaging (octaves)')
   ylabel('Correlation between ILDAlone & ILD/Freq Prediction');
   hold on
   for cell = 1:size(test_numbers,1)
      plot(cc_vec(cell,4),cc_bp_vec(cell),'b*')
   end
   xlabel('Correlation between ILDAlone & ILD/Freq Prediction, Tonal Avg')
   ylabel('Correlation between ILDAlone & ILD/Freq Prediction, BP Measured');
   
   for fig = 2:2:26
      figure(fig)
      print
   end

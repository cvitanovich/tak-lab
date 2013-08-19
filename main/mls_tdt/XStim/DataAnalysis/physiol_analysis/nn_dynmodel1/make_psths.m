%Script to produce psth's for the appropriate tests and save them to the Neuron structure

clear

load d:\mlspezio\matlab\save\Neuron_24_psth
clear ITD_matrix ILD_matrix ABI_matrix test_numbers


begin_cell = 1;
end_cell = 24;
start_time = 0.008; %seconds
end_time = .108; %seconds
Fs = 100; %Hz
bins = start_time:1/Fs:end_time-(1/Fs);

for cell_num = begin_cell:end_cell
   
   Neuron{cell_num}.psth_Fs = Fs;
   Neuron{cell_num}.psth_bins = bins;
   
   %Tonal psth
   clear data_array
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   data_array = Neuron{cell_num}.tif_dataarray;
   data_index = find(data_array(:,1) >= (start_time*1000) & data_array(:,1) <= (end_time*1000));
   clear data_array
   data_array = Neuron{cell_num}.tif_dataarray(data_index,:);
   data_array(:,1) = data_array(:,1)/1000;
   
   
   for freq_num = 1:length(Neuron{cell_num}.tif_freqaxis)
      for ild_num = 1:length(Neuron{cell_num}.tif_ildaxis)
         temp_index = find(data_array(:,4) == freq_num & data_array(:,5) == ild_num);
         [Neuron{cell_num}.tif_psth{freq_num,ild_num}] = hist(data_array(temp_index,1),bins);
      end
   end
	end
   
   %BP psth
   clear data_array
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
   data_array = Neuron{cell_num}.bif_dataarray;
   data_index = find(data_array(:,1) >= (start_time*1000) & data_array(:,1) <= (end_time*1000));
   clear data_array
   data_array = Neuron{cell_num}.bif_dataarray(data_index,:);
   data_array(:,1) = data_array(:,1)/1000;
   
   
   for freq_num = 1:length(Neuron{cell_num}.bif_freqaxis)
      for ild_num = 1:length(Neuron{cell_num}.bif_ildaxis)
         temp_index = find(data_array(:,4) == freq_num & data_array(:,5) == ild_num);
         [Neuron{cell_num}.bif_psth{freq_num,ild_num}] = hist(data_array(temp_index,1),bins);
      end
   end
	end
   
   %ILDAlone psth
   clear data_array
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ia_meansurf')));
   data_array = Neuron{cell_num}.ia_dataarray;
   data_index = find(data_array(:,1) >= (start_time*1000) & data_array(:,1) <= (end_time*1000));
   clear data_array
   data_array = Neuron{cell_num}.ia_dataarray(data_index,:);
   data_array(:,1) = data_array(:,1)/1000;
   
   
   for loc_num = 1:length(Neuron{cell_num}.ia_locs{1})
      ele = Neuron{cell_num}.ia_locs{1}(loc_num,1);
      azi = Neuron{cell_num}.ia_locs{1}(loc_num,2);
      temp_index = find(data_array(:,4) == ele & data_array(:,5) == azi);
      [Neuron{cell_num}.ia_psth{loc_num}] = hist(data_array(temp_index,1),bins);
   end
	end
   
   disp(['Finished processing cell # ' num2str(cell_num)])
end

save d:\mlspezio\matlab\save\Neuron_24_psth Neuron HRTFinfo
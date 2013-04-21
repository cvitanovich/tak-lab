function [isi,mean_isi,std_isi,Cv] = get_isi(data,time_col,rep_col)
%Function to process spike data and return InterSpikeInterval information
%data:					m X n, where n >= 2 and m are spikes
%time_col:				column in which time of spikes is recorded
%rep_col:				column in which repetition # of stimulus is recorded
%mean_isi:				mean ISI
%std_isi:				standard deviation of ISI
%Cv:						Coefficient of Variance (std/mean)

max_reps = max(data(:,rep_col));
diff_reps = diff(data(:,rep_col));
diff_reps = [-999;diff_reps];

zero_index = find(diff_reps == 0); %Only take those with diff==0 in order to group by rep & stim type
if(~isempty(zero_index))
   %for ind_num = 1:length(zero_index)
   %   isi(ind_num) = data(zero_index(ind_num),1) - data(zero_index(ind_num)-1,1);
   %end
   isi = data(zero_index,1) - data(zero_index-1,1);
   mean_isi = mean(isi);
   std_isi  = std(isi);
   Cv = std_isi/mean_isi;
else
   isi = 0;
   mean_isi = 0;
   std_isi = 0;
   Cv = 0;
end

return
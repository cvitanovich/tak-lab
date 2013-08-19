%Script to analyze the principal components of the ILD-distance relationship and
%determine if they correlate to a cell's tonal freq response
%	Fukunga, Introduction to Statistical Pattern Recognition (New York, NY: Academic Press, 1972)

clear
%set parameters
bird_number = 899;
cmap = jet;
linecolor = 'white';
bkgndcolor = [0.85 0.85 0.85];
%Directories & files
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\HRTF_files\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;
nntwarn off

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix

%Limit frequencies used
minfreq = 2000; maxfreq = 11000;
[y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
[y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
ILD_matrix = ILD_matrix(minfreqind:maxfreqind,:);

%Get those locations in the ILDmatrix that match the measured ILDAlone locations
for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
   ILDmat_index(num_loc) =...
      max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
      HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
end
ILD_matrix_focus = ILD_matrix(:,ILDmat_index);

for cell_num = 1:1
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ia_meansurf')));
      
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
      
      disp(['Processing cell # ' num2str(cell_num)])      
      
      %Calculate the distance matrix
      disp('Calculating ILD distances...')
      [dist_ild(:,:,cell_num)] = get_ilddist_all(ILD_matrix_focus,IA_bestloc_ind(cell_num));
      [dist_ild_signed(:,:,cell_num)] = get_ilddist_all_sign(ILD_matrix_focus,IA_bestloc_ind(cell_num));
      
      %Make the appropriate frequency-dependent activity data matrix
      disp('Making the activity data matrix...')
      [y,index] = sort(dist_ild_signed(:,:,cell_num),2);
      for freq_num = 1:size(dist_ild_signed,1)
         activity_mat(freq_num,:) = Neuron{cell_num}.ia_meansurf{1}(index(freq_num,:));
      end
      
      %Get Covariance Matrix
      disp('Calculating covariance matrix...')
      covmat = cov(activity_mat');
      
      %Get Principal Components
      disp('Finding principal components...')
      [princomp, latent, vari_explained] = pcacov(covmat);
      
      %Do the Karhunen-Loeve expansion
      disp('Calculating K-L transformation...')
      y_activity_mat(:,:,cell_num) = princomp'*...
         (activity_mat - repmat(mean(activity_mat,2),1,size(activity_mat,2)));
      
   end %end if
end %end loop over cells

disp('Finished all cells.')

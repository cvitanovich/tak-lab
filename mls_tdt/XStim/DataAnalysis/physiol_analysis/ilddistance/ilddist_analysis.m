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
%Define Frequency Bandlimits
freq_incr = 500;
fcount = 0;
for freq = 2000:freq_incr:10500
   fcount = fcount + 1;
   freq_bands(fcount,:) = [freq freq+freq_incr];
end

%Get those locations in the ILDmatrix that match the measured ILDAlone locations
for num_loc = 1:size(Neuron{1}.ia_locs{1},1)
   ILDmat_index(num_loc) =...
      max(find(HRTFinfo.location_matrix(1,:) == Neuron{1}.ia_locs{1}(num_loc,1) &...
      HRTFinfo.location_matrix(2,:) == Neuron{1}.ia_locs{1}(num_loc,2)));
end
ILD_matrix_focus = ILD_matrix(:,ILDmat_index);

for cell_num = 1:24
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
      
      
      [dist_ild(:,:,cell_num)] = get_ilddist(ILD_matrix_focus,...
         HRTFinfo.hrtf_freqs,...
         IA_bestloc_ind(cell_num),...
         freq_bands,0);
      
      fig_ht = 9.5; fig_wd = 7.5; mainfig = figure; set(mainfig,'Units','inches'); 
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      
      numcols = 3; numrows = round(size(freq_bands,1)/numcols);
      for fband = 1:size(freq_bands,1)
         h = subplot(numrows,numcols,fband);
         set(h,'FontSize',6);
         hold on
         plot(dist_ild(fband,:,cell_num),Neuron{cell_num}.ia_meansurf{1},'.');
         xlabel('Distance')
         ylabel('(spikes/s)')
         title(['(' num2str(freq_bands(fband,1))...
               ',' num2str(freq_bands(fband,2)) ')']);
      end
      
   end %end if
end %end loop over cells

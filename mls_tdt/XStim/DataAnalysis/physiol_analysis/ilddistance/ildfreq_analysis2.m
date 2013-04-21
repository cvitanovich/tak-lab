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
distflag = 1;

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix

for cell_num = 1:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
      
      %Find the best ILDAlone RS location
      IA_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ia_meansurf{1} == ...
         max(Neuron{cell_num}.ia_meansurf{1})));
      %if(cell_num == 1)  IA_bestloc_ind(cell_num) = 129; end
      %if(cell_num == 4)  IA_bestloc_ind(cell_num) = 198; end
      %if(cell_num == 7)  IA_bestloc_ind(cell_num) = 163; end
      %if(cell_num == 9)  IA_bestloc_ind(cell_num) = 282; end
      %if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
      %if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
      %if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
      %if(cell_num == 18) IA_bestloc_ind(cell_num) =  78; end
      %if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
      
      disp(['Processing cell # ' num2str(cell_num)])      
      clear ILD_matrix_focus ILD_matrix_new
      
		%Limit frequencies used
      minfreq = min(Neuron{cell_num}.tif_freqaxis);
      maxfreq = max(Neuron{cell_num}.tif_freqaxis);
      minILD = min(Neuron{cell_num}.tif_ildaxis);
      maxILD = max(Neuron{cell_num}.tif_ildaxis);
		[y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
		[y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
		ILD_matrix_focus = ILD_matrix(minfreqind:maxfreqind,:);
		freqs = HRTFinfo.hrtf_freqs(minfreqind:maxfreqind);
      
      %Get those locations in the ILDmatrix that match the measured ILDAlone locations
      for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
         ILDmat_index(num_loc) =...
            max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
            HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
      end
      ILD_matrix_focus = ILD_matrix_focus(:,ILDmat_index);
      
      %Restrict ILD matrix to only those frequencies in the Tonal measurement
      numfreqs = length(Neuron{cell_num}.tif_freqaxis);
      minfreq = min(Neuron{cell_num}.tif_freqaxis); maxfreq = max(Neuron{cell_num}.tif_freqaxis);
      freq_count = 0;
      for freq = minfreq:(maxfreq-minfreq)/(numfreqs-1):maxfreq
         freq_count = freq_count + 1;
         [y,freqind] = min(abs(freqs - freq));
         ILD_matrix_new(freq_count,:) = ILD_matrix_focus(freqind,:);
      end
      

      
%Get the frequency & ILD specifications to search for locations that share freq-spec ILD with
%IA_bestloc
      %Calculate the distance matrix
      disp('Calculating ILD distances...')
      [dist_ild] = get_ilddist_all(ILD_matrix_focus,IA_bestloc_ind(cell_num));
      clear cc_actfreq cc_actfreq_restrloc
      disp('Calculating correlations over all locations...')
      for n = 1:length(freqs)
      	temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},dist_ild(n,:));
         cc_actfreq(n) = temp(1,2);
      end
      %cc_actfreq = abs(cc_actfreq);
      cc_actfreq = -1*(cc_actfreq);
      [y,max_cc_ind] = max(abs(cc_actfreq));
      
      %Get the locations close to the optimal ILD_distance at the max_cc_ind
      clear ind_restrloc ind_high_act ind_low_act
      thresh_ILD_dist = 2; %+/- range in dB
      ind_restrloc = find( abs(ILD_matrix_focus(max_cc_ind,:) - ...
         ILD_matrix_focus(max_cc_ind,IA_bestloc_ind(cell_num))) <= thresh_ILD_dist);
      
      thresh_high_act = .8; thresh_low_act = .5;
      ind_high_act = find(Neuron{cell_num}.ia_meansurf{1}(ind_restrloc) >= ...
         thresh_high_act*max(Neuron{cell_num}.ia_meansurf{1}));
      ind_low_act  = find(Neuron{cell_num}.ia_meansurf{1}(ind_restrloc) <= ...
         thresh_low_act*max(Neuron{cell_num}.ia_meansurf{1}));
      
      disp('Calculating correlations over restricted locations...')
      for n = 1:length(freqs)
         temp = corrcoef(Neuron{cell_num}.ia_meansurf{1}(ind_restrloc),dist_ild(n,ind_restrloc));
         cc_actfreq_restrloc(n) = temp(1,2);
      end
      %cc_actfreq_restrloc = abs(cc_actfreq_restrloc);
      cc_actfreq_restrloc = -1*(cc_actfreq_restrloc);
      [y,max_cc_ind_restrloc] = max(abs(cc_actfreq_restrloc));
      
      
      
      %Plot
      fig_ht = 9.5; fig_wd = 7.5;      
      mainfig = figure; set(mainfig,'Units','inches'); 
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      subplot(3,1,1)
      plot(freqs,cc_actfreq,'k','LineWidth',2);
      set(gca,'FontSize',8);
      hold on
      plot(freqs,cc_actfreq_restrloc,'b-.');
      plot(Neuron{cell_num}.tif_freqaxis,...
         mean(Neuron{cell_num}.tif_meansurf',1)/max(mean(Neuron{cell_num}.tif_meansurf',1)),'m--');
      xlabel('Frequency (Hz)'); ylabel('Correlation (r)');
      title(['Cell # ' num2str(cell_num) ', with ' num2str(length(ind_restrloc))...
            ' locations +/- ' num2str(thresh_ILD_dist) ' dB at f = ' num2str(freqs(max_cc_ind))...
            ' Hz']);
      
      subplot(3,1,2)
      hist(Neuron{cell_num}.ia_meansurf{1}(ind_restrloc));
      
      
      %Measured Tonal ILD/freq RS
      subplot(3,1,3)
      set(gca,'FontSize',8); hold on;
      hp1 = plotsurf(Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         Neuron{cell_num}.tif_meansurf');
      ax1 = gca;
      xlim([minfreq maxfreq]); ylim([minILD maxILD]);
      colormap(cmap); colorbar;
      hold on
      if(exist('ind_high_act') & ~isempty(ind_high_act))
         plot(freqs,...
            ILD_matrix_focus(:,ind_restrloc(ind_high_act)),'g-');
      end
      if(exist('ind_low_act') & ~isempty(ind_low_act))
         plot(freqs,...
            ILD_matrix_focus(:,ind_restrloc(ind_low_act)),'m-.');
      end
      plot(freqs,...
         ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'LineWidth',2,'Color','white');
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      ax2 = axes('Position',get(ax1,'Position'));
      axes(ax2)
      hp2 = plot(freqs,cc_actfreq,'k--','LineWidth',2);
      set(ax2,...
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Color','none',...
         'XColor','k','YColor','k','FontSize',8);
      xlim([minfreq maxfreq]); ylim([0 0.8]); ylabel('Correlation (r)','Rotation',270);
      title(['Measured Tonal ILD/freq for cell # ' num2str(cell_num)]);
      
      
      
   end %end if
end %end loop over cells

disp('Finished all cells.')

%Script to analyze the principal components of the ILD-distance relationship and
%determine if they correlate to a cell's tonal freq response
%	Fukunga, Introduction to Statistical Pattern Recognition (New York, NY: Academic Press, 1972)

clear
%set parameters
bird_number = 899;
cmap = 1-gray;
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
      
      %Get ILD bins for Activity/ILD curve
      numILDs  = length(Neuron{cell_num}.tif_ildaxis);
      minILD  = min(Neuron{cell_num}.tif_ildaxis);  maxILD  = max(Neuron{cell_num}.tif_ildaxis);
      ILD_matrix_new(find(ILD_matrix_new > maxILD))  =  maxILD;
      ILD_matrix_new(find(ILD_matrix_new < minILD))  =  minILD;
      ILD_bins = minILD:(maxILD-minILD)/(numILDs-1):maxILD;

      
      
%Get the frequency & ILD specifications to search for locations that share freq-spec ILD with
%IA_bestloc
      %Calculate the distance matrix
      disp('Calculating ILD distances...')
      [dist_ild] = get_ilddist_all(ILD_matrix_focus,IA_bestloc_ind(cell_num));
      clear cc_actfreq;
      disp('Calculating correlations...')
      for n = 1:length(freqs)
      	temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},dist_ild(n,:));
         cc_actfreq(n) = temp(1,2);
      end
      cc_actfreq = abs(cc_actfreq);
      [y,max_cc_ind] = max(cc_actfreq);
      [y,freq_targ_ind] = min(abs(freqs(max_cc_ind) - Neuron{cell_num}.tif_freqaxis));
      [y,ild_targ_ind] = min(abs(ILD_bins - ILD_matrix_focus(max_cc_ind,IA_bestloc_ind(cell_num))));
      
      disp('Forming Activity/ILD curves for frequencies/locations')
      %Put the Activity/ILD relationship into ILD bins
      clear Activity; max_val = 0;
      for freq_num = 1:size(ILD_matrix_new,1)
         tempvec = ILD_matrix_new(freq_num,:);
         [n,ild_bin_index] = histc(tempvec,ILD_bins);
         for ild_bin_num = 1:length(ILD_bins)
            Activity(ild_bin_num,freq_num) = ...
               mean(Neuron{cell_num}.ia_meansurf{1}(find(ild_bin_index == ild_bin_num)));
            if(freq_num == freq_targ_ind & ild_bin_num == ild_targ_ind)
               temp_loc_ind = find(ild_bin_index == ild_bin_num);
            end
         end
      end
      
      count_ind = 0;
      count_ind2 = 0;
      thresh1 = 0.4; thresh2 = 0.7;
      clear key_loc_ind key_loc_ind2
      for n = 1:length(temp_loc_ind)
         if(Neuron{cell_num}.ia_meansurf{1}(temp_loc_ind(n)) < ...
               thresh1*max(Neuron{cell_num}.ia_meansurf{1}))
            count_ind = count_ind + 1;
            key_loc_ind(count_ind) = temp_loc_ind(n);
         end
         if(Neuron{cell_num}.ia_meansurf{1}(temp_loc_ind(n)) >= ...
               thresh2*max(Neuron{cell_num}.ia_meansurf{1}))
            count_ind2 = count_ind2 + 1;
            key_loc_ind2(count_ind2) = temp_loc_ind(n);
         end
      end
         
      
      %Change all NaN's to 0
      nanind = find(isnan(Activity) == 1);
      Activity(nanind) = 0;
      
      %Plot
      fig_ht = 9.5; fig_wd = 7.5;      
      mainfig = figure; set(mainfig,'Units','inches'); 
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      
      %Measured Tonal ILD/freq RS
      subplot(2,1,1)
      set(gca,'FontSize',8); hold on;
      hp1 = plotsurf(Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         Neuron{cell_num}.tif_meansurf');
      ax1 = gca;
      xlim([minfreq maxfreq]); ylim([minILD maxILD]);
      colormap(cmap); colorbar;
      hold on
      if(exist('key_loc_ind'))
         plot(freqs,...
            ILD_matrix_focus(:,key_loc_ind),'m-.');
      end
      if(exist('key_loc_ind2'))
         plot(freqs,...
            ILD_matrix_focus(:,key_loc_ind2),'g');
      end
      plot(freqs,...
         ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'LineWidth',2,'Color','blue');
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
      
      
      %Measured Narrow-Band Noise ILD/freq RS
      subplot(2,1,2)
      set(gca,'FontSize',8); hold on;
      if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
         plotsurf(Neuron{cell_num}.bif_freqaxis,Neuron{cell_num}.bif_ildaxis,...
            Neuron{cell_num}.bif_meansurf');
         colormap(cmap); colorbar;
         axis([minfreq maxfreq minILD maxILD]);
         hold on
      if(exist('key_loc_ind'))
         plot(freqs,...
            ILD_matrix_focus(:,key_loc_ind),'m-.');
      end
      if(exist('key_loc_ind2'))
         plot(freqs,...
            ILD_matrix_focus(:,key_loc_ind2),'g');
      end
         plot(freqs,...
            ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'LineWidth',2,'Color','blue');
      end
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      title('Measured Narrowband ILD/freq')
      
      
      
   end %end if
end %end loop over cells

disp('Finished all cells.')

%Analyze distance using local regression to fit Activity/ILD plots

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
distflag = 0;

%Local Regression Parameters
alpha = 0.05; mode = 2;

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

for cell_num = 1:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ia_meansurf')));
      
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
      
      %Calculate the distance matrix
      disp('Calculating ILD distances & sorting...')
      [dist_ild(:,:,cell_num)] = get_ilddist_all(ILD_matrix_focus,IA_bestloc_ind(cell_num));
      [dist_ild_sort(:,:,cell_num),dsort_index] = sort(dist_ild(:,:,cell_num),2);
      [dist_ild_signed(:,:,cell_num)] = get_ilddist_all_sign(ILD_matrix_focus,IA_bestloc_ind(cell_num));
      [dist_ild_s_sort(:,:,cell_num),dsort_s_index] = sort(dist_ild_signed(:,:,cell_num),2);
      
      [ILD_matrix_focus_sort,ILD_sort_index] = sort(ILD_matrix_focus,2);
      
      disp('Calculating correlations...')
      freqs = HRTFinfo.hrtf_freqs(minfreqind:maxfreqind);
      for n = 1:length(freqs)
      	temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},dist_ild(n,:,cell_num));
      	cc_actfreq(n,cell_num) = temp(1,2);
      end
      
      %Define the tonal activity and its frequency axis
      if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
         tone_ridge = max(Neuron{cell_num}.tif_meansurf',[],1);
         tif_freqaxis = Neuron{cell_num}.tif_freqaxis;
      else
         tif_freqaxis = zeros(1,21);
         tone_ridge = zeros(1,21);
      end
      
      %Make all correlations positive so that comparison between elements is easier
      cc_actfreq = abs(cc_actfreq);
      
      %Get frequency index of maximal correlation between Cell Activity and ILD distance
      [y,max_ccafi_ind] = max(cc_actfreq(:,cell_num));
      
      %Get frequency indices surrounding the index of maximal correlation
      freq_range = 4000; %Hz
      for n = 1:5
         switch n
         case {1,2},
            [y,freq_ind_vec(n)] = min(abs(freqs - (freqs(max_ccafi_ind) - freq_range/(2*n))));
         case 3,
            freq_ind_vec(n) = max_ccafi_ind;
         case {4,5},
            [y,freq_ind_vec(n)] = min(abs(freqs - (freqs(max_ccafi_ind) + freq_range/(2*(6-n)))));
         end
      end
      
      %Use local regression (locregress1) to approximate Activity vs. (ILD or ILD_dist) function
      disp('Finding best fits by local regression...')
      if(distflag) %for ILD *distance*
         for n = 1:5
            [result(n,:),residual(n,:)] = locregress1(dist_ild_s_sort(freq_ind_vec(n),:,cell_num),...
               Neuron{cell_num}.ia_meansurf{1}(dsort_s_index(freq_ind_vec(n),:)),alpha,mode);
            %Interpolate to bring the function closer to the peak
         end
      else %for plain ILD
         for n = 1:5
            [result(n,:),residual(n,:)] = locregress1(ILD_matrix_focus_sort(freq_ind_vec(n),:),...
               Neuron{cell_num}.ia_meansurf{1}(ILD_sort_index(freq_ind_vec(n),:)),alpha,mode);
         end
      end
         
      
      %Plot all
      disp('Plotting...')
      
      %Plot the Activity vs. ILD distance for frequencies specified by max corr freq and range
      fig_ht = 9.5; fig_wd = 7.5;      
      secondfig = figure; set(secondfig,'Units','inches'); 
      set(secondfig,'Position',[0 0 fig_wd fig_ht]);
      count = 0;
      hand = subplot(2,3,1);
      set(gca,'FontSize',8); hold on;
      [haxes,hline1,hline2] = plotyy(freqs,cc_actfreq(:,cell_num),tif_freqaxis,tone_ridge,'plot');
      set(hline1,'LineStyle','-','Color','k','LineWidth',2);
      set(hline2,'LineStyle','--','Color','m','LineWidth',2);
      axes(haxes(1)); xlabel('Frequency (Hz)');
      set(gca,'YColor','black');
      axes(haxes(1)); ylabel('Correlation (r)');
      line(freqs(freq_ind_vec(1))*ones(1,10),...
         0:max(cc_actfreq(:,cell_num))/9:max(cc_actfreq(:,cell_num)),'LineWidth',1,'Color','black');
      line(freqs(freq_ind_vec(5))*ones(1,10),...
         0:max(cc_actfreq(:,cell_num))/9:max(cc_actfreq(:,cell_num)),'LineWidth',1,'Color','black');
      axes(haxes(2)); ylabel('Activity (spikes/s)','FontSize',7,'Rotation',270);
      set(gca,'FontSize',8,'YColor','magenta');
      title(['r (solid), Activity (dashed), cell # ' num2str(cell_num)]);
      for n = 1:5
         freq_ind = freq_ind_vec(n);
         count = count + 1;
         figure(secondfig);
         hand = subplot(2,3,count+1);
         set(hand,'FontSize',8)
         if(distflag)
            plot(dist_ild_signed(freq_ind,:,cell_num),Neuron{cell_num}.ia_meansurf{1},'k.');
            hold on
            plot(dist_ild_s_sort(freq_ind,:,cell_num),result(n,:),'LineWidth',1.5);
            xlabel('Distance from ILDAlone optimum (dB)');
         else
            plot(ILD_matrix_focus(freq_ind,:),Neuron{cell_num}.ia_meansurf{1},'k.');
            hold on
            plot(ILD_matrix_focus_sort(freq_ind,:),result(n,:),'LineWidth',1.5);
            xlabel('ILD (dB)');
         end
         ylabel('Activity (spikes/stimulus)');
         axis([-20 20 min(Neuron{cell_num}.ia_meansurf{1}) max(Neuron{cell_num}.ia_meansurf{1})]);
         title(['Cell # ' num2str(cell_num) ' @ ' num2str(round(freqs(freq_ind))) ' Hz'])
      end
         
      
   end %end if
end %end loop over cells

disp('Finished all cells.')

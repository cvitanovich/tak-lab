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

fig_ht = 9.5; fig_wd = 7.5; mainfig = figure; set(mainfig,'Units','inches'); 
set(mainfig,'Position',[0 0 fig_wd fig_ht]);
plotcount = 0;
for cell_num = 22:22
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
            
      %Get Covariance Matrix
      disp('Calculating covariance matrix...')
      covmat = cov(dist_ild(:,:,cell_num)');
      
      %Get Principal Components
      disp('Finding principal components...')
      [princomp, latent, vari_explained] = pcacov(covmat);
      
      disp('Calculating K-L transformation...')
      y_dist_ild(:,:,cell_num) = princomp'*...
         (dist_ild(:,:,cell_num) - repmat(mean(dist_ild(:,:,cell_num),2),1,size(dist_ild(:,:,cell_num),2)));
      
      disp('Calculating correlations...')
      freqs = HRTFinfo.hrtf_freqs(minfreqind:maxfreqind);
      for n = 1:length(freqs)
         temp = corrcoef(dist_ild(n,:,cell_num),y_dist_ild(1,:,cell_num));
         cc_ilddist(n,1,cell_num) = temp(1,2);
         temp = corrcoef(dist_ild(n,:,cell_num),y_dist_ild(2,:,cell_num));
         cc_ilddist(n,2,cell_num) = temp(1,2);
         temp = corrcoef(dist_ild(n,:,cell_num),y_dist_ild(3,:,cell_num));
         cc_ilddist(n,3,cell_num) = temp(1,2);
         
      	temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},dist_ild(n,:,cell_num));
      	cc_actfreq(n,cell_num) = temp(1,2);
      end
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},y_dist_ild(1,:,cell_num));
      cc_actpc(1,cell_num) = temp(1,2);
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},y_dist_ild(2,:,cell_num));
      cc_actpc(2,cell_num) = temp(1,2);
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},y_dist_ild(3,:,cell_num)); 
      cc_actpc(3,cell_num) = temp(1,2);
      
      %Define the tonal activity and its frequency axis
      if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
         tone_ridge = max(Neuron{cell_num}.tif_meansurf',[],1);
         tif_freqaxis = Neuron{cell_num}.tif_freqaxis;
      else
         tif_freqaxis = zeros(1,21);
         tone_ridge = zeros(1,21);
      end
      
      %Make all correlations positive so that comparison between elements is easier
      cc_ilddist = abs(cc_ilddist);
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
      
      %Minimize the corresponding ILD_distance function
      disp('Minimizing Activity/ILD_distance function')
   	Options = optimset('Display','iter','MaxIter',100,'TolFun',0.002,'TolX',0.01,'MaxFunEvals',1000);
      V = [17 0 1/17 0]; %width of curve, height of sigmoidal portion
      [V,err,exitflag] = ...
         fminsearch('sig_ilddistfun',...
         V,...
         Options,...
         Neuron{cell_num}.ia_meansurf{1},...
         dist_ild_signed(max_ccafi_ind,:,cell_num));
      
      [err,result,V] = ...
         sig_ilddistfun(V,...
         Neuron{cell_num}.ia_meansurf{1},...
         dist_ild_signed(max_ccafi_ind,:,cell_num));
      [y,plotindex] = sort(dist_ild_signed(max_ccafi_ind,:,cell_num));
      result = result(plotindex);
      ildwidth(cell_num) = V(1); err_all_lin(cell_num) = err;
      
      %Plot all
      disp('Plotting...')
      
      clear haxes hline1 hline2
      numplots = 6;
      plotcount = plotcount + 1;
      figure(mainfig);
      subplot(3,2,plotcount)
      set(gca,'FontSize',8); hold on;
      hold on;
      [haxes,hline1,hline2] = plotyy(freqs,cc_ilddist(:,1,cell_num),...
         tif_freqaxis,tone_ridge,'plot');
      set(hline1,'LineStyle','-','Color','b','LineWidth',1);
      set(hline2,'LineStyle','--','Color','m','LineWidth',2);
      axes(haxes(1)); xlabel('Frequency (Hz)');
      set(gca,'YColor','black');
      axes(haxes(1)); ylabel('Correlation (r)');
      axes(haxes(2)); ylabel('Activity (spikes/s)','FontSize',7,'Rotation',270);
      set(gca,'FontSize',8,'YColor','magenta');
      axes(haxes(1));
      plot(freqs,cc_ilddist(:,2,cell_num),'r-.','LineWidth',2);
      plot(freqs,cc_ilddist(:,3,cell_num),'g:','LineWidth',2);
      plot(freqs,cc_actfreq(:,cell_num),'k','LineWidth',3);
      axes(haxes(2));
      if(0)
         legend(['ILDdist/1PC, ' num2str(vari_explained(1)) ' % var'],...
            ['ILDdist/2PC, ' num2str(vari_explained(2)) ' % var'],...
            ['ILDdist/3PC, ' num2str(vari_explained(3)) ' % var'],...
            'Act/ILDdist',...
            'Norm Tonal Max',-1);
      end
      title(['Cell # ' num2str(cell_num)]);
      
      if(mod(cell_num,numplots) == 0 & cell_num < length(Neuron))
         fig_ht = 9.5; fig_wd = 7.5; mainfig = figure; set(mainfig,'Units','inches'); 
         set(mainfig,'Position',[0 0 fig_wd fig_ht]);
         plotcount = 0;
      end
      
      %Plot the Activity vs. ILD distance for frequencies specified by max corr freq and range
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
         plot(dist_ild_signed(freq_ind,:,cell_num),Neuron{cell_num}.ia_meansurf{1},'k.');
         if(n == 3)
            hold on;
            plot(sort(dist_ild_signed(freq_ind,:,cell_num)),result,'r-','LineWidth',2);
         end
         axis([-40 40 min(Neuron{cell_num}.ia_meansurf{1}) max(Neuron{cell_num}.ia_meansurf{1})]);
         if(n == 3)
            xlabel('Distance from ILDAlone Optimum')
            ylabel('Activity (spikes/s)')
         end
         title(['Cell # ' num2str(cell_num) ' @ ' num2str(round(freqs(freq_ind))) ' Hz'])
      end
         
      
   end %end if
end %end loop over cells

disp('Finished all cells.')

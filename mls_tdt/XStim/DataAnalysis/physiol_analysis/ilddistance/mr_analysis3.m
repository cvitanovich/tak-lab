%mr_analysis3
%Compare tonal and bb multiple regression results, using the weights determined by 
%multiple regression in both cases

clear;

%set parameters
bird_number = 899;
begin_cells = 1;
tot_num_cells = 24;
end_cells = begin_cells + tot_num_cells - 1;

%Graphics parameters
bkgndcolor = [0.85 0.85 0.85];
fig_ht = 9.5;
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
      ILDf_azi			= Neuron{cell_num}.tia_azi;
      ILDf_ele			= Neuron{cell_num}.tia_ele;
      ILDf_diamond	= Neuron{cell_num}.tia_diamond;
      
      IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
      IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
      IA_diamond  = Neuron{cell_num}.ia_diamond{1};
      IA_locs		= Neuron{cell_num}.ia_locs{1};
      IA_azi		= Neuron{cell_num}.ia_azi{1};
      IA_ele		= Neuron{cell_num}.ia_ele{1};
      
      if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')))
         TS_meansurf = Neuron{cell_num}.ts_meansurf{1};
         TS_azi      = Neuron{cell_num}.ts_azi{1};
         TS_ele      = Neuron{cell_num}.ts_ele{1};
         TS_diamond  = Neuron{cell_num}.ts_diamond{1};
      end
      
      %Get tonal ridge
      ILDf_ridge = max(ILDf_meansurf,[],2);
      
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
   
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')))
      TS_bestloc_ind(cell_num) = max(find(TS_meansurf == max(TS_meansurf)));
   end
   
   %Load up the multiple regression results
   %1. Tonal
   eval(['load d:\mlspezio\matlab\save\cell' num2str(cell_num) '_mrstats']);
   tonal_in = in; tonal_beta = zeros(size(beta));
   tonal_beta(in) = beta(in); clear in beta;
   %2. BB
   eval(['load d:\mlspezio\matlab\save\cell' num2str(cell_num) '_mrstats1a']);
   bb_in = in; bb_beta = zeros(size(beta));
   bb_beta(in) = beta(in); clear in beta;
   
   %Calculate tonal ILD-alone RS
   tia_meansurf = 0;
   for freq_num = tonal_in'
      temp_ILDf_meansurf = zeros(size(ILDf_meansurf));
      temp_ILDf_meansurf(freq_num,:) = ILDf_meansurf(freq_num,:);
      [temp,dirs,full_rs,new_f_ax,new_ild_ax] = ...
         if2space(ILD_matrix_focus',...
         HRTFinfo.hrtf_freqs,...
         ILDmat_index,...
         temp_ILDf_meansurf,...
         ILDf_freqaxis,...
         ILDf_ildaxis);
      %This step makes sure that we don't have completely 0 rows
      if(~any(temp))
         temp(1) = 0.01;
      end
      tia_meansurf = tia_meansurf + (tonal_beta(freq_num) * temp);
      disp(['Finished tonal simulation for ' num2str(ILDf_freqaxis(freq_num)) ' Hz'])
   end
   cc_mrtone = corrcoef(IA_meansurf,tia_meansurf);
   var_mrtone = cc_mrtone(1,2)^2;
   
   %Calculate the bb ILD-alone RS
   %Determine Frequency Bands to use in the MR analysis
   for freq_num = 1:length(ILDf_freqaxis)-1
      freq_bands(freq_num,1) = ILDf_freqaxis(freq_num);
      freq_bands(freq_num,2) = ILDf_freqaxis(freq_num + 1);
   end
   %Calculate distance_ILD matrix, unsigned
   disp('Calculating signed distance matrix...')
   sign_flag = 1;
   dist_ild_signed = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
      IA_bestloc_ind(cell_num),freq_bands,sign_flag);
   [dist_ild_s_sort,dsort_s_index] = sort(dist_ild_signed,2);
   %***Important: reverse the sorting for the comparison
   [junk,reverse_sort_index] = sort(dsort_s_index,2);
   clear junk
   %Use local regression (locregress1) to approximate Activity vs. (ILD or ILD_dist) function
   disp('Finding best fits by local regression...')
   alpha = 0.15;
   clear tempia resia bbia_meansurf
   bbia_meansurf = 0;
   for freq_num = bb_in'
      [tempia,resia] = locregress1(dist_ild_s_sort(freq_num,:),...
         IA_meansurf(dsort_s_index(freq_num,:)),alpha);
      bbia_meansurf = bbia_meansurf + (bb_beta(freq_num)* tempia(reverse_sort_index(freq_num,:)));
      disp(['Finished bb simulation for ' num2str(ILDf_freqaxis(freq_num)) ' Hz'])
      %subplot(10,2,freq_num)
      %plot(dist_ild_s_sort(freq_num,:),tempia);
   end
   cc_bb = corrcoef(IA_meansurf,bbia_meansurf);
   var_bb = cc_bb(1,2)^2;

   %Plot
   h = figure;
   set(h,'Units','inches','Position',[0 0 fig_wd fig_ht]);
   
   %Tonal ridge
   subplot('Position',[0.05 0.8 0.45 0.16]);
   plot(ILDf_freqaxis,ILDf_ridge)
   hold on
   plot(ILDf_freqaxis(tonal_in'),ILDf_ridge(tonal_in'),'b*')
   plot([ILDf_freqaxis(bb_in') ILDf_freqaxis(bb_in'+1)]',...
      [ILDf_ridge(bb_in') ILDf_ridge(bb_in'+1)]','r','Linewidth',2)
   xlabel('Frequency (Hz)'); ylabel('Activity (sp/stim)');
   title(['Cell # ' num2str(cell_num) ', ridge of ILD/frequency RS']);
   axis([min(ILDf_freqaxis) max(ILDf_freqaxis) 0 max(ILDf_ridge)]);
   str1(1) = {'*  : Tonal MR'}; str1(2) = {'---: BB MR'};
   text(max(ILDf_freqaxis),max(ILDf_ridge),str1(1),'HorizontalAlignment','Left','Color','b');
   text(max(ILDf_freqaxis),max(ILDf_ridge)-2,str1(2),'HorizontalAlignment','Left','Color','r');
   
   %ILD/frequency surface
   subplot('Position',[0.61 0.8 0.37 0.16]);
   plotsurf(ILDf_freqaxis,ILDf_ildaxis,ILDf_meansurf');
   hold
   plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'y','Linewidth',1.5);
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')))
      plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,TS_bestloc_ind(cell_num)),'w--','Linewidth',1.5);
   end
   xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
   title('ILD/frequency RS')
   
   %Space
   subplot('Position',[0.07 0.50 0.2 0.2])
   axis square
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')))
      plotdiam(TS_azi,TS_ele,TS_diamond);
   end
   xlabel('Azimuth (deg)'); ylabel('Elevation (deg)'); title('Spatial RS');
   
   %ILD-alone
   subplot('Position',[0.31 0.50 0.2 0.2])
   plotdiam(IA_azi,IA_ele,IA_diamond);
   xlabel('Azimuth (deg)'); title('ILD-alone RS');
   
   %Tonal ILD-alone, Linear
   subplot('Position',[0.55 0.50 0.2 0.2])
   plotdiam(ILDf_azi,ILDf_ele,ILDf_diamond);
   xlabel('Azimuth (deg)');
   cc_lintone = corrcoef(IA_meansurf,ILDf_meanarray(ILDmat_index));
   var_lintone = cc_lintone(1,2)^2;
   title(str2mat('Tonal ILD-alone RS',['LINEAR, r^{2} = ' num2str(var_lintone)]),'FontSize',8);
   
   %Tonal ILD-alone, MR
   subplot('Position',[0.79 0.50 0.2 0.2])
   %make the diamond for the prediction
   [IA_mrtonepredict_azi,IA_mrtonepredict_ele,temp] = ...
      array2diamond(tia_meansurf,IA_locs');
	%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(IA_mrtonepredict_azi, IA_mrtonepredict_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(IA_locs(:,1)==el & IA_locs(:,2)==az))
	       missmask(AZ==az & EL==el) = 1;
	    else
	       missmask(AZ==az & EL==el) = 0;
	    end; 
	  end;
	end;
	
	% replace all valid locations presently containing NaN's with zeros
	ri = (missmask==1);
	temp(ri) = zeros(size(temp(ri)));
	% generate interpolation function - uses surrounding 4 squares to generate missing square
   intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
	temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
	intval_old = conv2(temp, intfun,'same')/4;
	intval = intval_old.*missmask;
   IA_mrtonepredict_diam = intval + temp;
   plotdiam(IA_mrtonepredict_azi,IA_mrtonepredict_ele,IA_mrtonepredict_diam);
   xlabel('Azimuth (deg)');
   title(str2mat('Tonal ILD-alone RS',['MR, r^{2} = ' num2str(var_mrtone)]),'FontSize',8);
   
   %BB ILD-alone, MR
   subplot('Position',[0.07 0.25 0.2 0.2])
   %make the diamond for the prediction
   [IA_mrbbpredict_azi,IA_mrbbpredict_ele,temp] = ...
      array2diamond(bbia_meansurf,IA_locs');
	%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(IA_mrbbpredict_azi, IA_mrbbpredict_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(IA_locs(:,1)==el & IA_locs(:,2)==az))
	       missmask(AZ==az & EL==el) = 1;
	    else
	       missmask(AZ==az & EL==el) = 0;
	    end; 
	  end;
	end;
	
	% replace all valid locations presently containing NaN's with zeros
	ri = (missmask==1);
	temp(ri) = zeros(size(temp(ri)));
	% generate interpolation function - uses surrounding 4 squares to generate missing square
   intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
	temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
	intval_old = conv2(temp, intfun,'same')/4;
	intval = intval_old.*missmask;
   IA_mrbbpredict_diam = intval + temp;
   plotdiam(IA_mrbbpredict_azi,IA_mrbbpredict_ele,IA_mrbbpredict_diam);
   xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
   title(str2mat('ILD_{dist}  ILD-alone RS',['MR, r^{2} = ' num2str(var_bb)]),'FontSize',8);
   
end %end if
end %end for
   
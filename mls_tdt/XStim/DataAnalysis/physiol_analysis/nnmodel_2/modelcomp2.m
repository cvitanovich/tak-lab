%Script to compare Linear distance-based model with nnmodel_2

clear;
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

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix
load d:\mlspezio\matlab\save\vars_1_24;
load d:\mlspezio\matlab\save\dbp_elements;

for cell_num = 1:1
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
      
		neg_index = find(bp_optimal_ildf_surface{cell_num} < 0);
		bp_optimal_ildf_noinhib{cell_num} = bp_optimal_ildf_surface{cell_num};
		bp_optimal_ildf_noinhib{cell_num}(neg_index) = 0;
      numfreqs = length(Neuron{cell_num}.tif_freqaxis);
      numILDs  = length(Neuron{cell_num}.tif_ildaxis);
      
      %Contour for the spatial receptive field of the neuron
      if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
         [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
            srfstats(Neuron{cell_num}.ts_azi{1},...
            Neuron{cell_num}.ts_ele{1},...
            Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
      end
      
      % eliminate location of AZ=0, EL=90 so that it won't give us problems
      oklocs = HRTFinfo.location_matrix(1,:)~=90;
      HRTFinfo.location_matrix = HRTFinfo.location_matrix(:,oklocs);
      ILD_matrix = ILD_matrix(:,oklocs);
      
      %****************************
      %Linmodel9
      %****************************
      %Get those locations in the ILDmatrix that match the measured ILDAlone locations
      for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
         ILDmat_index(num_loc) =...
            max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
            HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
      end
      ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
      
      %Get the frequency indices
      if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
         [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - Neuron{cell_num}.tif_freqaxis(1)));
         [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
            Neuron{cell_num}.tif_freqaxis(length(Neuron{cell_num}.tif_freqaxis))));
      end
      
      %Find the best ILDAlone RS location
      IA_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ia_meansurf{1} == max(Neuron{cell_num}.ia_meansurf{1})));
      if(cell_num == 1) IA_bestloc_ind(cell_num) = 129; end
      if(cell_num == 4) IA_bestloc_ind(cell_num) = 198; end
      if(cell_num == 7) IA_bestloc_ind(cell_num) = 163; end
      if(cell_num == 9) IA_bestloc_ind(cell_num) = 282; end
      if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
      if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
      if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
      if(cell_num == 18) IA_bestloc_ind(cell_num) = 78; end
      if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
      
      %Interpolate the ILD/Freq surface at the resolution of the HRTF-derived ILD spectra
      [new_if_surface,new_if_freqaxis,ILDmat_freqind] = ...
         if_interp(Neuron{cell_num}.tif_meansurf,...
         Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         HRTFinfo.hrtf_freqs);
      
      %Define Frequency Bandlimits
      freq_incr = 1000;
      fcount = 0;
      for freq = 2000:freq_incr:10000
         fcount = fcount + 1;
         freq_bands(fcount,:) = [freq freq+freq_incr];
      end
      if(isempty(freq_bands))
         size_freq_bands = size(new_if_freqaxis,1);
      else
         size_freq_bands = size(freq_bands,1);
      end
      
      %Calculate Model prediction with bestfit values V
      if(cell_num < 21)
         row_ind = cell_num;
      else
         row_ind = cell_num - 1;
      end
      [err,result,V,freq_weight] = errfun_specint1(vars_1_24(row_ind,:),...
         Neuron{cell_num}.ia_meansurf{1},...
         new_if_surface,...
         new_if_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         ILD_matrix_focus(ILDmat_freqind(1):ILDmat_freqind(2),:),...
         IA_bestloc_ind(cell_num),...
         freq_bands);
      
      IA_linpredict = result;
      index_neg = find(IA_linpredict < 0);
      IA_linpredict(index_neg) = 0;
      
      %make the diamond for the prediction
      [IA_linpredict_azi,IA_linpredict_ele,temp] = ...
         array2diamond(IA_linpredict,Neuron{cell_num}.ia_locs{1}');
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(IA_linpredict_azi, IA_linpredict_ele);
      
      % generate mask for missing points
      missmask = NaN*ones(size(temp));
      
      i = 1;
      for az = -90:5:90;
         for el = -90+abs(az):5:90-abs(az)
            if (~(Neuron{cell_num}.ia_locs{1}(:,1)==el & Neuron{cell_num}.ia_locs{1}(:,2)==az))
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
      
      IA_linpredict_diam = intval + temp;
      
      %Calculate Predicted ILDAlone RS from Euston Backprop - NO INHIBITION
      [IA_dbpmodel_ni_meansurf, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
         HRTFinfo.hrtf_freqs,...
         HRTFinfo.location_matrix,...
         bp_optimal_ildf_noinhib{cell_num},...
         Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis);
      [IA_dbpmodel_ni_azi,IA_dbpmodel_ni_ele,IA_dbpmodel_ni_diamond] = ...
         array2diamond(IA_dbpmodel_ni_meansurf,HRTFinfo.location_matrix);
      
      %Calculate Predicted Tonal ILD/Freq surface using BroadBand minimized values
      [TonalILDf_pred,V,freq_weight] = ...
         mk_ildfreq(V,...
         Neuron{cell_num}.tif_meansurf,...
         Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         ILD_matrix_focus,...
         HRTFinfo.hrtf_freqs,...
         IA_bestloc_ind(cell_num),...
         freq_bands);
      
      clear temp
      temp = corrcoef(Neuron{cell_num}.tif_meansurf(:),TonalILDf_pred(:));
      Tonal_cc(cell_num) = temp(1,2)^2;
      
      %PLOTTING
      fig_ht = 9.5;
      fig_wd = 7.5;
      mainfig = figure;
      set(mainfig,'Units','inches');
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      
      %Measured ILDAlone RS
   	subplot('Position',[0.07 0.65 0.25 0.25])
      set(gca,'FontSize',8);
      axis square
      plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},Neuron{cell_num}.ia_diamond{1});
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      title(['Measured ILDAlone for Cell #' num2str(cell_num)],'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Predicted TONAL ILDAlone RS
   	subplot('Position',[0.07 0.35 0.25 0.25])
      set(gca,'FontSize',8);
      axis square   
      plotdiam(Neuron{cell_num}.tia_azi,Neuron{cell_num}.tia_ele,Neuron{cell_num}.tia_diamond);
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      %Find the correlation coefficient between the Tonal ILD/freq predicted and measured ILDAlone RS
      for cc_loc = 1:length(Neuron{cell_num}.ia_locs{1})
         cc_index(cc_loc) = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,1)...
            & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,2)));
      end
      if(any(strcmp(fieldnames(Neuron{cell_num}),'tia_meanarray')));
         temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},Neuron{cell_num}.tia_meanarray(cc_index));
         cc = temp(1,2)^2;
         ccTonal(cell_num) = cc;
      end
      title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc)]),'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Predicted David_Backprop_Model ILDAlone RS
      subplot('Position',[0.38 0.07 0.25 0.25])
      set(gca,'FontSize',8);
      axis square
      plotdiam(IA_dbpmodel_azi,IA_dbpmodel_ele,IA_dbpmodel_diam(:,:,cell_num));
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},squeeze(IA_dbpmodel_meansurf(1,ILDmat_index,cell_num)));
      cc = temp(1,2)^2;
      ccDBPmodel(cell_num) = cc;
      title(str2mat('Euston Backprop Model',['r^{2} = ', num2str(cc)]),'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Predicted David_Backprop_Model ILDAlone RS - NO INHIBITION
      subplot('Position',[0.71 0.07 0.25 0.25])
      set(gca,'FontSize',8);
      axis square
      plotdiam(IA_dbpmodel_ni_azi,IA_dbpmodel_ni_ele,IA_dbpmodel_ni_diamond);
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},IA_dbpmodel_ni_meansurf(ILDmat_index));
      cc = temp(1,2)^2;
      ccDBPmodel_ni(cell_num) = cc;
      title(str2mat('Euston Backprop Model - NO INHIBITION',['r^{2} = ', num2str(cc)]),'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Predicted Linmodel9 ILDAlone RS
   	subplot('Position',[0.07 0.07 0.25 0.25])
      set(gca,'FontSize',8);
      axis square
      plotdiam(IA_linpredict_azi,IA_linpredict_ele,IA_linpredict_diam);
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},IA_linpredict);
      cc = temp(1,2)^2;
      ccLinmodel(cell_num) = cc;
      title(str2mat('Linear-model Predicted ILDAlone',['r^{2} = ', num2str(cc)]),'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Select frequency range matching the tif_meansurf
      minfreq = min(Neuron{cell_num}.tif_freqaxis); maxfreq = max(Neuron{cell_num}.tif_freqaxis);
      [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
      [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
      
      ILD_matrix_new = interp1(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),...
         ILD_matrix_focus(ind_minfreq:ind_maxfreq,:),...
         Neuron{cell_num}.tif_freqaxis,'spline');
      
      %Measured ILD/freq RS
   	subplot('Position',[0.4 0.7 0.55 0.23])
      set(gca,'FontSize',8);
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
      hold on
      colorbar
      colormap(cmap);
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'yellow',...
         'LineWidth',2);
      axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
            min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
      xlabel('Freq (Hz)','FontSize',8);
      ylabel('ILD (dB)','FontSize',8);
      title('Measured Tonal ILD/Freq','FontSize',8);
      
      %Predicted ILD/freq RS
   	subplot('Position',[0.4 0.4 0.55 0.23])
      set(gca,'FontSize',8);
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,TonalILDf_pred);
      hold on
      colorbar
      colormap(cmap);
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'yellow',...
         'LineWidth',2);
      axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
            min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
      xlabel('Freq (Hz)','FontSize',8);
      ylabel('ILD (dB)','FontSize',8);
      title(['Predicted Tonal ILD/Freq, r^{2} = ' num2str(Tonal_cc(cell_num))],'FontSize',8);
      
      
      
      
      
      clear P T space_contour
      
   end %end asking if Neuron has a TIF
end %end loop over neurons

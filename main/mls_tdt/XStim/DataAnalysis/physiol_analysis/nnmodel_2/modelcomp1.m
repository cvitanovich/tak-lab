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
load d:\mlspezio\matlab\save\nnet_24b;
%nnet_24{1}.IW{1}(25,25) = 5; nnet_24{1}.IW{1}(305,305) = 2;nnet_24{1}.IW{1}(311,311) = 2;
%nnet_24{1}.IW{1}(212,212) = 2;

for cell_num = 24:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
      
      numfreqs = length(Neuron{cell_num}.tif_freqaxis);
      numILDs  = length(Neuron{cell_num}.tif_ildaxis);
      
      %Contour for the spatial receptive field of the neuron
      if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
         [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
            srfstats(Neuron{cell_num}.ts_azi{1},...
            Neuron{cell_num}.ts_ele{1},...
            Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
      end
      
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
      [err,result,V,freq_weight] = errfun_specint1(vars_1_24(23,:),...
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
      
      %****************************
      %nnmodel_2
      %****************************
      
      %Select frequency range matching the tif_meansurf
      minfreq = min(Neuron{cell_num}.tif_freqaxis); maxfreq = max(Neuron{cell_num}.tif_freqaxis);
      [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
      [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
      
      ILD_matrix_new = interp1(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),...
         ILD_matrix_focus(ind_minfreq:ind_maxfreq,:),...
         Neuron{cell_num}.tif_freqaxis,'spline');
      
      %Set input ranges
      minILD = min(min(ILD_matrix_new)); maxILD = max(max(ILD_matrix_new));
      ICspacenet{cell_num}.inputs{1}.range = [minILD*ones(numfreqs*numILDs,1) maxILD*ones(numfreqs*numILDs,1)];
      
      %Set inputs and targets - each input/target pair is an ILD spectrum for a location
      %and the unit's response to that location
      [FREQs,ILDs] = meshgrid(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis);
      tempP = zeros(length(Neuron{cell_num}.tif_ildaxis),length(Neuron{cell_num}.tif_freqaxis),size(ILD_matrix_new,2));
      for loc = 1:size(ILD_matrix_new,2)
         for freq = 1:size(ILD_matrix_new,1)
            [y,ildindex]  = min(abs(Neuron{cell_num}.tif_ildaxis  - ILD_matrix_new(freq,loc)));
            tempP(ildindex,freq,loc) = 1;
         end
         %P is the input we use. Each column is a location containing numfreqs*numILDs values
         %1's are in the appropriate matrix cells for a location's ILD/freq spectrum
         P(:,loc) = reshape(tempP(:,:,loc),...
            length(Neuron{cell_num}.tif_ildaxis)*length(Neuron{cell_num}.tif_freqaxis),1);
      end
      TonalP = zeros(numfreqs*numILDs,numfreqs*numILDs);
      for inputnum = 1:size(TonalP,2)
         TonalP(inputnum,inputnum) = 1;
      end
      P = [P TonalP]; %Include Tonal inputs for training set
      T = Neuron{cell_num}.ia_meansurf{1}; %ILDAlone targets
      b = ones(1,numfreqs*numILDs);
      b(:) = Neuron{cell_num}.tif_meansurf';
      T = [T b]; clear b;
      IA_nnpredict = sim(nnet_24b{cell_num},P(:,1:359));
      [IA_nnpredict_azi,IA_nnpredict_ele,temp] = ...
         array2diamond(IA_nnpredict,Neuron{cell_num}.ia_locs{1}');
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(IA_nnpredict_azi, IA_nnpredict_ele);
      
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
      IA_nnpredict_diam = intval + temp;
      
      %PLOTTING
      fig_ht = 9.5;
      fig_wd = 7.5;
      mainfig = figure;
      set(mainfig,'Units','inches');
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      
      %Measured ILDAlone RS
      subplot(3,3,1)
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
      
      %Measured ILD/freq RS
      subplot(3,3,3)
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
      
      %Predicted TONAL ILDAlone RS
      subplot(3,3,2)
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
         cctonal(cell_num) = cc;
      end
      title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc)]),'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Predicted Linmodel9 ILDAlone RS
      subplot(3,3,4)
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
      
      %Predicted nnmodel_2 ILDAlone RS
      subplot(3,3,5)
      set(gca,'FontSize',8);
      axis square
      plotdiam(IA_nnpredict_azi,IA_nnpredict_ele,IA_nnpredict_diam);
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},IA_nnpredict);
      cc = temp(1,2)^2;
      cc_nnmodel(cell_num) = cc;
      title(str2mat('NN-model Predicted ILDAlone',['r^{2} = ', num2str(cc)]),'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Predicted nnmodel_2 ILD/freq RS
      subplot(3,3,6)
      set(gca,'FontSize',8);
      z = sim(nnet_24b{cell_num},TonalP);
      negind = find(z < 0); z(negind) = 0;
      tif_nnpredict = zeros(numILDs,numfreqs);
      tif_nnpredict(:) = z;
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,tif_nnpredict);
      hold on
      colorbar
      colormap(cmap);
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'yellow',...
         'LineWidth',2);
      axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
            min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
      xlabel('Freq (Hz)','FontSize',8);
      ylabel('ILD (dB)','FontSize',8);
      title('Modeled Tonal ILD/Freq','FontSize',8);
      
      %Trained input weights
      subplot(3,3,7)
      set(gca,'FontSize',8);
      a = zeros(numILDs,numfreqs);
      a(:) = diag(nnet_24b{cell_num}.IW{1});
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,a);
      shading flat
      colorbar
      colormap(cmap)
      xlabel('Freq (Hz)','FontSize',8); ylabel('ILD (dB)','FontSize',8);
      title('Input weights, after training','FontSize',8)
      
      %Trained Layer 1-->2 weights
      subplot(3,3,8)
      set(gca,'FontSize',8);
      a = zeros(numILDs,numfreqs);
      a(:) = nnet_24b{cell_num}.LW{2,1};
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,a);
      shading flat
      colorbar
      colormap(cmap)
      xlabel('Freq (Hz)','FontSize',8); ylabel('ILD (dB)','FontSize',8);
      title('Weights for layer 1-->2, after training','FontSize',8)
      hold on
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'k',...
         'LineWidth',2);
      
      clear P T space_contour
   end %end asking if Neuron has a TIF
end %end loop over neurons

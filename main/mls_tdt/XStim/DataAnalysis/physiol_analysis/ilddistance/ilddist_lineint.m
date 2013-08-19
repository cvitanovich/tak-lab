%Script to test Line Integral method WITH ILD_distance approach

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

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'

for cell_num = 1:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
      
      disp(['Processing cell # ' num2str(cell_num)])
      
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
      
      %Get the contour for the SRF of the neuron 
      if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
         [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
            srfstats(Neuron{cell_num}.ts_azi{1},...
            Neuron{cell_num}.ts_ele{1},...
            Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
      end
      
      %Get those locations in the ILDmatrix that match the measured ILDAlone locations
      for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
         ILDmat_index(num_loc) =...
            max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
            HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
      end
      ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
      
      
      disp('Getting predicted ILDAlone RS with Line Integrat/ILDdist method...')
      %Get Predicted ILDAlone RS
      clear ILDf_mat_all;
      [ILDf_mat_all] = if2space3(ILD_matrix_focus,...
         HRTFinfo.hrtf_freqs,...
         Neuron{cell_num}.tif_meansurf',...
         Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis);
      
      disp('Transforming the matrix...')
      %Transform matrix into freq X loc matrix
      clear New_ILDf_mat
      for loc_num = 1:size(ILDf_mat_all,3)
         for freq_num = 1:length(Neuron{cell_num}.tif_freqaxis)
            ildind = find(ILDf_mat_all(:,freq_num,loc_num) ~= 0);
            if(isempty(ildind))
               New_ILDf_mat(freq_num,loc_num) = 0;
            else
               New_ILDf_mat(freq_num,loc_num) = ILDf_mat_all(ildind,freq_num,loc_num);
            end
         end
      end
      
      %minimize the model error
      disp('Minimizing parameters...')
   	Options = optimset('Display','iter','MaxIter',8000,'TolFun',0.002,'TolX',0.01,'MaxFunEvals',30000);
      
      clear V
      V = 10*ones(1,length(Neuron{cell_num}.tif_freqaxis));
      [V,err] = fminsearch('ilddist_weights',V,Options,...
         New_ILDf_mat,...
         ILD_matrix_focus,...
         HRTFinfo.hrtf_freqs,...
         Neuron{cell_num}.tif_meansurf',...
         Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         IA_bestloc_ind(cell_num),...
         Neuron{cell_num}.ia_meansurf{1});
      
      lineint_V{cell_num} = V;
      save d:\mlspezio\matlab\save\lineint_V lineint_V 
      clear IA_tonalpredict
      [err,IA_tonalpredict] = ilddist_weights(V,...
         New_ILDf_mat,...
         ILD_matrix_focus,...
         HRTFinfo.hrtf_freqs,...
         Neuron{cell_num}.tif_meansurf',...
         Neuron{cell_num}.tif_freqaxis,...
         Neuron{cell_num}.tif_ildaxis,...
         IA_bestloc_ind(cell_num),...
         Neuron{cell_num}.ia_meansurf{1});

      
      %make the diamond for the prediction
      [IA_tonalpredict_azi,IA_tonalpredict_ele,temp] = ...
         array2diamond(IA_tonalpredict,Neuron{cell_num}.ia_locs{1}');
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(IA_tonalpredict_azi, IA_tonalpredict_ele);
      
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
      
      IA_tonalpredict_diam = intval + temp;
      
      %Calculate the variance explained
      clear temp
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},IA_tonalpredict);
      vari_exp = temp(1,2)^2;
      
      %Plot
      fig_ht = 9.5;
      fig_wd = 7.5;
      mainfig = figure;
      set(mainfig,'Units','inches');
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      
      subplot('Position',[.07 .07 .55 .25]);
      set(gca,'FontSize',8);
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
      hold on
      plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),'w','LineWidth',2);
      axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis) ...
            min(Neuron{cell_num}.tif_ildaxis) max(Neuron{cell_num}.tif_ildaxis)]);
      xlabel('Frequency (Hz)');
      ylabel('ILD (dB)');
      title('Measured ILD/freq RS');
      
      subplot('Position',[.07 .4 .25 .25]);
      set(gca,'FontSize',8);
      axis square
      plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},Neuron{cell_num}.ia_diamond{1});
      colormap(cmap); colorbar;
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),'white','LineWidth',1.5);
      end
      xlabel('Azimuth');
      ylabel('Elevation');
      title(['Measured ILDAlone RS, cell # ' num2str(cell_num)]);
      
      subplot('Position',[.4 .4 .25 .25]);
      set(gca,'FontSize',8);
      axis square
      plotdiam(Neuron{cell_num}.tia_azi,Neuron{cell_num}.tia_ele,Neuron{cell_num}.tia_diamond);
      colormap(cmap); colorbar;
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),'white','LineWidth',1.5);
      end
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},Neuron{cell_num}.tia_meanarray(ILDmat_index));
      vari_exp_tonal(cell_num) = temp(1,2)^2;
      xlabel('Azimuth');
      ylabel('Elevation');
      title(str2mat('Tonal Line Integral ILDAlone RS',['r^{2} = ' num2str(vari_exp)]));
      
      subplot('Position',[.7 .4 .25 .25]);
      set(gca,'FontSize',8);
      axis square
      plotdiam(IA_tonalpredict_azi,IA_tonalpredict_ele,IA_tonalpredict_diam);
      colormap(cmap); colorbar;
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),'white','LineWidth',1.5);
      end
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},IA_tonalpredict);
      vari_exp_tonal_ilddist(cell_num) = temp(1,2)^2;
      xlabel('Azimuth');
      ylabel('Elevation');
      title(str2mat('Tonal Line Integral/ILDdist ILDAlone RS',['r^{2} = ' num2str(vari_exp)]));
      
      save d:\mlspezio\matlab\save\vari_exp_000729 vari_exp_tonal vari_exp_tonal_ilddist
   end %end if
end %end for

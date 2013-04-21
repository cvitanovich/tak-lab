%Script to generate figures for Initial Planning of Poster for ARO2001 meeting



%Load neuronal data
if(~exist('Neuron','var'))
   load 'd:\mlspezio\matlab\save\Neuron_50_regress_twelfth'
end

%Load HRTF files for 901, 913 & 894
if(exist('d:\mlspezio\matlab\save\HRTFAll_4.mat'))
   eval(['load d:\mlspezio\matlab\save\HRTFAll_4'])
end

new_ildaxis = [-30:1:30];
vert_pos = [0.50 0.05];


for cell_num = [1]
   
   clear ILD_matrix HRTFinfo ILDf_freqaxis ILDf_ildaxis ILDf_meansurf ILDf_stdsurf ILD_matrix_focus
   clear IA_meansurf IA_stdsurf IA_diamond IA_locs IA_azi IA_ele bbif_meansurf beta pc pc_bbif_meansurf
   
   %Specify Site Parameters
   if(cell_num < 25)
      bird_number = 899;
      depth = Neuron{cell_num}.tif_testpars(7);
      site_number = cell_num;
      side_of_brain = '';
      stim_dur = Neuron{cell_num}.tif_testpars(8);
      isi = [];
      session_num = [];
   else
      bird_number = Neuron{cell_num}.abi_testpars.bird_number;
      session_num = Neuron{cell_num}.abi_testpars.session_num;
      depth = Neuron{cell_num}.abi_testpars.depth;
      site_number = Neuron{cell_num}.abi_testpars.recording_site;
      side_of_brain = Neuron{cell_num}.abi_testpars.side_of_brain;
      stim_dur = Neuron{cell_num}.abi_testpars.curr_stimdur;
      isi = Neuron{cell_num}.abi_testpars.test_ISI;
   end
   
   %Choose HRTF file
   if(cell_num < 25) %for 899
      if(0)
         hrtf_file = 'd:\mlspezio\owl_hrtfdata\899\out9be';
         get_HRTF = 1;
         get_ITDmatrix = 0;
         hrtf_file2 = [];
         [ILD_matrix,ABL_matrix,ITD_matrix,HRTFinfo] = ...
            get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
         save 'd:\mlspezio\matlab\save\899HRTF' ILD_matrix ABL_matrix HRTFinfo
      else
         load 'd:\mlspezio\matlab\save\899HRTF' ILD_matrix ABL_matrix HRTFinfo
      end
   else
      switch bird_number
      case 901
         ILD_matrix = HRTFAll{1}.ILD_matrix;
         HRTFinfo = HRTFAll{1}.HRTFinfo;
      case 913
         ILD_matrix = HRTFAll{2}.ILD_matrix;
         HRTFinfo = HRTFAll{2}.HRTFinfo;
      case 894
         ILD_matrix = HRTFAll{3}.ILD_matrix;
         HRTFinfo = HRTFAll{3}.HRTFinfo;
      end %end switch
   end %end if
   
   %Calculate ILD-distance-based measures
      disp(['Processing cell # ' num2str(cell_num)])
      
      %Choose the elements for analysis:
      if(isfield(Neuron{cell_num},'gif_freqaxis'))
         ILDf_freqaxis  = Neuron{cell_num}.gif_freqaxis;
         ILDf_ildaxis   = Neuron{cell_num}.gif_ildaxis;
         ILDf_meansurf  = Neuron{cell_num}.gif_mean;
         ILDf_meanarray = Neuron{cell_num}.gia_meanarray;
      elseif(isfield(Neuron{cell_num},'tif_freqaxis'))
         ILDf_freqaxis  = Neuron{cell_num}.tif_freqaxis;
         ILDf_ildaxis   = Neuron{cell_num}.tif_ildaxis;
         ILDf_meansurf  = Neuron{cell_num}.tif_meansurf;
         ILDf_meanarray = Neuron{cell_num}.tia_meanarray;
      end
      
      if(cell_num < 25)
         IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
         IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
         IA_diamond  = Neuron{cell_num}.ia_diamond{1};
         IA_locs		= Neuron{cell_num}.ia_locs{1};
         IA_azi		= Neuron{cell_num}.ia_azi{1};
         IA_ele		= Neuron{cell_num}.ia_ele{1};
      else
         IA_meansurf = Neuron{cell_num}.ia_mean{1};
         IA_stdsurf  = Neuron{cell_num}.ia_std{1};
         IA_diamond  = Neuron{cell_num}.ia_diamond{1};
         IA_locs		= Neuron{cell_num}.ia_locs{1};
         IA_azi		= Neuron{cell_num}.ia_azi{1};
         IA_ele		= Neuron{cell_num}.ia_ele{1};
      end
      
      if(size(IA_locs,1) < 3)
         IA_locs = IA_locs';
      end
      
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
      
      if(isfield(Neuron{cell_num},'ts_meansurf'))
         TS_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ts_meansurf{1} == max(Neuron{cell_num}.ts_meansurf{1})));
      elseif(isfield(Neuron{cell_num},'ts_mean'))
         TS_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ts_mean{1} == max(Neuron{cell_num}.ts_mean{1})));
      end
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
      
      
      %*****Begin Analysis for 1/12-octave bands
      
      %Determine Frequency Bands to use in the MR analysis
      clear freq_bands
      freq_bands = get_freqbands([2000 11000],1/12);
      bb_ILDf_freqaxis = freq_bands(:,1);
      
      %Calculate distance_ILD matrix, unsigned and signed
      disp('Calculating signed distance matrix...')
      sign_flag = 1;
      dist_ild_signed = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
         IA_bestloc_ind(cell_num),freq_bands,sign_flag);
      dist_ild{cell_num} = dist_ild_signed;
      [dist_ild_s_sort,dsort_s_index] = sort(dist_ild_signed,2);
      %***Important: reverse the sorting for the comparison
      [junk,reverse_sort_index] = sort(dsort_s_index,2);
      clear junk
      
      %Calculate Activity vs. ILD_distance curves
      alpha = 0.15;
      numpts = 10;
      correction_flag = 0;
      %Correct for wild tails if necessary
      %corr_cell_nums = [2 7 9 10 13 18 19 23 25 39];
      %if(~all(corr_cell_nums - cell_num)) correction_flag = 1;end
      thresh = 0.3;
      for freq_num = 1:size(freq_bands,1)
         [tempia,resia] = locregress1(dist_ild_s_sort(freq_num,:),...
            IA_meansurf(dsort_s_index(freq_num,:)),alpha,0);
         %[tempia,resia] = Actd_est(dist_ild_s_sort(freq_num,:),...
         %   IA_meansurf(dsort_s_index(freq_num,:)),numpts);
         len = length(tempia);
         %correct for tails that go wild
         maxval = 0.5*max(tempia);
         beg_ind = find(tempia(1:3) > maxval);
         end_ind = find(tempia(len-2:len) > maxval);
         if(~isempty(beg_ind) & correction_flag) %Correct tail at beginning
            [y,beg_min_ind] = min(tempia(1:50));
            if(beg_min_ind > 0)
               tempia(1:beg_min_ind) = interp1([dist_ild_s_sort(freq_num,1) dist_ild_s_sort(freq_num,beg_min_ind)],...
                  [min(tempia) tempia(beg_min_ind)],dist_ild_s_sort(freq_num,1:beg_min_ind));
            end
         end
         if(~isempty(end_ind) & correction_flag) %Correct tail at end
            [y,end_min_ind] = min(tempia(len-49:len));
            end_min_ind = len - (50 - end_min_ind);
            if(end_min_ind < len-3)
               tempia(end_min_ind:len) = interp1([dist_ild_s_sort(freq_num,end_min_ind) dist_ild_s_sort(freq_num,len)],...
                  [tempia(end_min_ind) min(tempia)],dist_ild_s_sort(freq_num,end_min_ind:len));
            end
         end
         tempia_allfreqs(freq_num,:) = tempia;
         resia_allfreqs(freq_num,:) = resia;
         bbif_meansurf(freq_num,:) = tempia(reverse_sort_index(freq_num,:));
         disp(['Finished bb simulation for ' num2str(freq_bands(freq_num,1)) '-' num2str(freq_bands(freq_num,2)) ' Hz'])
      end
      
      %Perform PCR (Principal Components Regression) on bbif_meansurf
      [pc,score,latent,tsq] = princomp(bbif_meansurf');
      frac_latent = latent/sum(latent);
      temp = find(frac_latent >= 0.05) %Keep all principal components contributing >= threshhold variance
      npc = length(temp) %number of PCs kept
      pc_totvar(cell_num) = sum(frac_latent(temp));
      
      %Transform bbif_meansurf by its PCs
      for pc_num = 1:npc
         pc_bbif_meansurf(:,pc_num) = (sum(repmat(pc(:,pc_num),1,size(bbif_meansurf,2)) .* bbif_meansurf,1))';
      end
      
      %Perform multiple regression using the PCs
      alpha = 0.05;
      if(0)
         stepwise(pc_bbif_meansurf,IA_meansurf',1:size(pc_bbif_meansurf,2),alpha);
         keyboard
      else
         [b1, delta1, tstats1, stats1] = regress_mls(pc_bbif_meansurf,IA_meansurf',[1:size(pc_bbif_meansurf,2)],alpha);        
         inmodel = [];
         bint = [b1-delta1 b1+delta1];
         for b_num = 1:length(b1)
            if(sign(bint(b_num,1)) == sign(bint(b_num,2))) inmodel(b_num) = b_num; end
         end
         inmodel = dezero(inmodel);
      end
      [b2, delta, tstats, stats] = regress_mls(pc_bbif_meansurf,IA_meansurf',inmodel,alpha);
      beta = zeros(size(b2));
      beta(inmodel) = b2(inmodel);
      
      final_coef = pc(:,1:npc) * beta; %scale the PC's with the regression coefficients
      
      Test_meanarray = final_coef' * bbif_meansurf;
      temp = corrcoef(IA_meansurf,Test_meanarray);
      pervar(cell_num) = temp(1,2)^2

      
      %Transform to scaled coefficients
      Y = IA_meansurf; X = bbif_meansurf;
      syy = sum((Y - mean(Y)).^2);
      for xnum = 1:size(bbif_meansurf,1)
         sxx(xnum) = sum((X(xnum,:) - mean(X(xnum,:))).^2);
         scale_coef(xnum) = final_coef(xnum) * sqrt(sxx(xnum)/syy);
      end
      
      norm_scale_coef = scale_coef/max(abs(scale_coef));
      
      %Interpolate to fit Activity/ILD curves onto newildaxis
      if(size(bb_ILDf_freqaxis,2) > 1) bb_ILDf_freqaxis = bb_ILDf_freqaxis'; end
      new_bbif_meansurf{cell_num} = griddata(repmat(bb_ILDf_freqaxis,1,size(bbif_meansurf,2)),...
         dist_ild_signed,...
         bbif_meansurf,...
         bb_ILDf_freqaxis,...
         new_ildaxis,'nearest');
      
      new_bbif_meansurf{cell_num} = (new_bbif_meansurf{cell_num})';
      
      %Make sure we're not interpolating out of bounds
      for freq_num = 1:length(bb_ILDf_freqaxis)
         tempind1 = find(new_ildaxis < min(dist_ild_signed(freq_num,:)));
         tempind2 = find(new_ildaxis > max(dist_ild_signed(freq_num,:)));
         new_bbif_meansurf{cell_num}(freq_num,[tempind1 tempind2]) = NaN;
      end
      
      %Weight by the Normalized, Scaled Coefficients
      for freq_num = 1:length(bb_ILDf_freqaxis)
         new_bbif_meansurf{cell_num}(freq_num,:) = ...
            new_bbif_meansurf{cell_num}(freq_num,:) * norm_scale_coef(freq_num);
      end
      
      %Transform by correcting for signed RMS ILD at optimal location
      [new_bbif_meansurf2{cell_num},new_ildaxis2] = ...
         iftransform_regress(Neuron,...
         (new_bbif_meansurf{cell_num}),...
         bb_ILDf_freqaxis,...
         new_ildaxis,...
         ILD_matrix_focus(:,IA_bestloc_ind(cell_num)),...
         HRTFinfo);
      
      [y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - min(ILDf_freqaxis)));
      [y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - max(ILDf_freqaxis)));
      
      hf = figure;
      set(hf,'Units','inches',...
         'Position',[0 0 7 5]);
      %List information for the Neuron
      Neuron_info = [num2str(bird_number) ...
            side_of_brain ...
         num2str(session_num) ','...
            ' Rec Site '...
            num2str(site_number) ','...
            ' Depth = '...
            num2str(depth) ','...
            ' Stim Dur = '...
            num2str(stim_dur) ','...
            ' ISI = '...
            num2str(isi)];
      
      htext = uicontrol('Parent',hf,...
         'Style','text',...
         'Units','normal',...
         'Position',[0 0.98 1 0.02],...
         'String',Neuron_info,...
         'FontSize',10,...
         'FontWeight','bold');
      
      %TIF or GIF
      hildf1 = subplot('Position',[0.05 vert_pos(1) 0.90 0.35]);
      set(hildf1,'FontSize',8);
      if(isfield(Neuron{cell_num},'gif_freqaxis'))
         plotsurf(Neuron{cell_num}.gif_freqaxis,Neuron{cell_num}.gif_ildaxis,Neuron{cell_num}.gif_mean');
         title('ILD-Freq, Gammatones');
      else
         plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
         title('ILD-Freq, Tones');
      end
      if(isfield(Neuron{cell_num},'ia_azi'))
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind(cell_num)),...
            'm',...
            'LineWidth',1.5);
      end
      if(isfield(Neuron{cell_num},'ts_azi'))
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind(cell_num)),...
            'w',...
            'LineWidth',1.5);
      end
      xlim([2000 11000]);
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      
      
      %BBIF
      hildf2 = subplot('Position',[0.05 vert_pos(2) 0.90 0.35]);
      set(hildf2,'FontSize',8);
      plotsurf(bb_ILDf_freqaxis,new_ildaxis2,(new_bbif_meansurf2{cell_num})');
      if(isfield(Neuron{cell_num},'ia_azi'))
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind(cell_num)),...
            'm',...
            'LineWidth',1.5);
      end
      if(isfield(Neuron{cell_num},'ts_azi'))
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind(cell_num)),...
            'w',...
            'LineWidth',1.5);
      end
      xlim([2000 11000]);
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      
   end %end for loop over cells
   
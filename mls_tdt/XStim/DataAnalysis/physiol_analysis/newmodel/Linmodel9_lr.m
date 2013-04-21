%Linmodel9 - script descibing and testing a linear model of spectral integration in ILDAlone RS
clear;
Model = 'Linmodel9_{LR}';

%set parameters
bird_number = 899;
begin_cells = 1;
end_cells = 24;
pca_weights = 0;
temporal_flag = 1;

%Graphics parameters
cmap = jet;
linecolor = 'white';
bkgndcolor = [0.85 0.85 0.85];
fig_ht = 10;
fig_wd = 8;
Dot_markersize = 4;
plot_dots = 0; %flag for plotting

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

%Local Regression Parameters
alpha_dist = 0.05; alpha_ILD = 0.05; mode = 0;

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'

%Test numbers
test_numbers = [
%[Thr  SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]
  -57  -1  -1  65    -1  66  68  70  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  71  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -62  -1  74  75    -1  76  77  79  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  80  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1  -1  85    -1  87  88  89  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  90  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1  -1  91    -1  92  93  94  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  95  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -80  -1  96  97    -1  98 100 101  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1 103 104    -1 105 106 107  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 108 109    -1 110 111 112  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 113  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 121 122    -1 123 124 125  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 126  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1  -1 127    -1 129 130 131  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 132  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 133 134    -1 135 136 137  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 138  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 141 142    -1 143 144 145  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 146  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 152 153    -1 162 158 159  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 161  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 165 166    -1 167 168 169  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 170  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 171 172    -1 173 174 175  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 176  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 206 208    -1 209 210 211  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 207 212 213  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70 220 214 216   223 217 218 219  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 221  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 228 230   234 232 233 231  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 229  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 242 248    -1 249 250 247  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 246  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 264 265    -1 274  -1 267 268  -1 269  -1 270 271 272 273   -1   -1   -1 266  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 289 290    -1 300  -1 293 294  -1 295  -1 296 297 298 299   -1   -1   -1 292  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -60  -1 318 319    -1  -1  -1 321 322  -1 323  -1  -1  -1  -1  -1   -1   -1   -1 320  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 326 327    -1 332  -1 329 330  -1 331 333 334 335  -1  -1   -1   -1   -1 328  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 337 338    -1 343  -1 340 341  -1 342  -1  -1  -1  -1  -1   -1   -1   -1 339  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 347 348    -1 356  -1 351 352  -1 353  -1  -1 355 358 357   -1   -1   -1 349  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
];

%get the HRTF info
%[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = ...
%   get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
%clear ITD_matrix ABI_matrix

%Begin looping over the neurons
for cell_num = begin_cells:end_cells
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
      
      %Choose the elements for analysis:
      %1. ILD_frequency frequency axis
      %2. ILDAlone meansurf and diamond
      
      if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
         ILDf_freqaxis = Neuron{cell_num}.bif_freqaxis;
      else
         ILDf_freqaxis = Neuron{cell_num}.tif_freqaxis;
      end
      
      if(temporal_flag == 1)
         time_limits = [8 60];
         [IA_meansurf,IA_stdsurf,IA_locs,IA_diamond,IA_azi,IA_ele] = ...
            get_IA_temporal(bird_number,side_of_brain,test_numbers(cell_num,8),1,0,time_limits);
      else
         IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
         IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
         IA_diamond  = Neuron{cell_num}.ia_diamond{1};
         IA_locs		= Neuron{cell_num}.ia_locs{1};
         IA_azi		= Neuron{cell_num}.ia_azi{1};
         IA_ele		= Neuron{cell_num}.ia_ele{1};
      end
      
   disp(['Processing cell # ' num2str(cell_num)])
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
   
   %Get the frequency indices
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      ILDf_freqaxis(length(ILDf_freqaxis))));
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
   
   
   disp(['Modelling activity for ' num2str(length(Neuron{cell_num}.ia_meansurf{1})) ' locations...'])
   clear dist_ild dist_ild_signed dist_ild_signed_sort d_s_sort_index result_pieces
   clear residual result ActILD_result ActILD_residual ILD_matrix_focus_restr_sort ILDf_linpredict
   clear IA_linpredict ILD_bands freq_bands
   bandwidth = 1/3; %units of octaves
   
   %Determine Frequency Bands (from Tonal ILD/frequency
   for freq_num = 1:length(ILDf_freqaxis)
      %Determine frequency limits
      [freq_bands(freq_num,1),freq_bands(freq_num,2)] = ...
         bandlimits(ILDf_freqaxis(freq_num),bandwidth);
   end
   
   %Calculate the RMS ILD across the frequency bands
   disp('Getting ILD across frequency bands...')
   sign_flag = 1;
   [ILD_bands] = get_ild_bands(ILD_matrix_focus,HRTFinfo.hrtf_freqs,freq_bands,sign_flag);

   %Calculate distance_ILD matrix, unsigned and signed
   disp('Calculating distance matrices...')
   sign_flag = 0;
   dist_ild = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
      IA_bestloc_ind(cell_num),freq_bands,sign_flag);
   sign_flag = 1;
   dist_ild_signed = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
      IA_bestloc_ind(cell_num),freq_bands,sign_flag);
   
   disp('Sorting ILD and distance matrices...')
   %Sort ILD matrix and distance matrices
   %[ILD_matrix_focus_sort,ILD_sort_index] = sort(ILD_matrix_focus,2);
   [ILD_bands_sort,ILD_bands_sort_index] = sort(ILD_bands,2);
   [dist_ild_sort,di_sort_index] = sort(dist_ild,2);
   [dist_ild_signed_sort,dis_sort_index] = sort(dist_ild_signed,2);
   
   if(plot_dots)
      firstfig = figure('Units','inches','Position',[0 0 fig_wd fig_ht],'PaperPositionMode','auto');
      numfreqs = length(ILDf_freqaxis);
      numfigcols = 3; numfigrows = ceil(numfreqs/numfigcols);
      
      secondfig = figure('PaperPositionMode','auto');
      set(gcf,'Units','inches');
      set(gcf,'Position',[0 0 fig_wd fig_ht]);
   end
   
   disp('Fitting Activity/ILD and Activity/distance curves...')
   for freq_num = 1:length(ILDf_freqaxis)
      [y,freqind] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(freq_num)));
      
      [ActILD_temp(freq_num,:),ActILD_residual{freq_num}] = ...
         locregress1(ILD_bands_sort(freq_num,:),...
         Neuron{cell_num}.ia_meansurf{1}(ILD_bands_sort_index(freq_num,:)),alpha_ILD,mode);
            
      %Calculate correlation between Activity & ILD_distance
      temp = corrcoef(dist_ild_sort(freq_num,:),Neuron{cell_num}.ia_meansurf{1}(di_sort_index(freq_num,:)));
      cc_actfreq(freq_num) = -1*temp(1,2);
      
      %Set Frequency weights
      cc_threshhold = 0;
      if(cc_actfreq(freq_num) >= cc_threshhold) 
         freq_weights(freq_num) = cc_actfreq(freq_num);
      else
         freq_weights(freq_num) = 0;
      end
      
      %Perform Local Regression to get the Activity vs. ILD_distance functions
      [Actilddist(freq_num,:),result_residual{freq_num}] = locregress1(dist_ild_signed_sort(freq_num,:),...
         Neuron{cell_num}.ia_meansurf{1}(dis_sort_index(freq_num,:)),alpha_dist,mode);
      
      [junk,reverse_sort_index(freq_num,:)] = sort(dis_sort_index,2);
      clear junk
      
      result_pieces(freq_num,:) = ...
         freq_weights(freq_num) * Actilddist(freq_num,reverse_sort_index(freq_num,:));
      ActILD_result(freq_num,:) = freq_weights(freq_num) * ...
         ActILD_temp(freq_num,:);
      
      if(plot_dots)
         figure(firstfig)
         subplot(numfigrows,numfigcols,freq_num)
         set(gca,'FontSize',6,'FontName','Garamond');
         plot(ILD_bands_sort(freq_num,:),...
            Neuron{cell_num}.ia_meansurf{1}(ILD_bands_sort_index(freq_num,:)),'k.',...
            'MarkerSize',Dot_markersize);
         hold on
         plot(ILD_bands_sort(freq_num,:),ActILD_temp(freq_num,:),'g','LineWidth',1.5);
         plot(Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf(freq_num,:),'r--','LineWidth',1.5);
         %plot(ILD_matrix_focus_restr_sort(freq_num,:),ActILD_temp,'b*','LineWidth',1.5);
         title(['Freq = ' num2str(ILDf_freqaxis(freq_num))],...
            'FontName','Garamond');
         xlim([min(Neuron{cell_num}.tif_ildaxis) max(Neuron{cell_num}.tif_ildaxis)])
         axes('Units','Normal','Position',[0 0 1 1]);
         text(0.5,0.98,['Cell # ' num2str(cell_num) ', Activity/ILD'],...
            'HorizontalAlignment','Center','FontName','Garamond','FontSize',12);
         set(gca,'Visible','off');
      
         figure(secondfig)
         subplot(numfigrows,numfigcols,freq_num)
         set(gca,'FontSize',6,'FontName','Garamond');
         plot(dist_ild_signed_sort(freq_num,:),...
            Neuron{cell_num}.ia_meansurf{1}(dis_sort_index(freq_num,:)),'k.',...
            'MarkerSize',Dot_markersize);
         hold on
         plot(dist_ild_signed_sort(freq_num,:),...
            Actilddist(freq_num,:),'m','LineWidth',1.5);
         %plot(ILD_matrix_focus_restr_sort(freq_num,:),ActILD_temp,'b*','LineWidth',1.5);
         title(['Freq = ' num2str(ILDf_freqaxis(freq_num))],...
            'FontName','Garamond');
         xlim([min(Neuron{cell_num}.tif_ildaxis) max(Neuron{cell_num}.tif_ildaxis)])
         axes('Units','Normal','Position',[0 0 1 1]);
         text(0.5,0.98,['Cell # ' num2str(cell_num) ', Activity/ILD_{dist}'],...
            'HorizontalAlignment','Center','FontName','Garamond','FontSize',12);
         set(gca,'Visible','off');
      end
   end
   
   
   %Choose between 2 ways of weighting the frequency bands
   if(pca_weights)
      %Find principal components of the Actilddist, treating ILD (i.e., location) as an observation
      disp('Calculating principal components...')
      [pcomp{cell_num},eigvals{cell_num},vari_exp{cell_num},kl_trans{cell_num}] = ...
         pca_kl(Actilddist');
      
      freq_weights2 = 0;
      for c_num = 1:3
         freq_weights2 = freq_weights2 + ((vari_exp{cell_num}(c_num)/100)*pcomp{cell_num}(:,c_num));
      end
      freq_weights2 = abs(freq_weights2);
      
      clear ActILD_result2
      for freq_num = 1:length(ILDf_freqaxis)
         result_pieces2(freq_num,:) = freq_weights2(freq_num) * ...
            Actilddist(freq_num,reverse_sort_index(freq_num,:));
         ActILD_result2(freq_num,:) = freq_weights2(freq_num) * ...
            ActILD_temp(freq_num,:);
      end
      result = sum(result_pieces2,1);
      
      %Calculate modelled ILDf meansurf
      disp('Calculating ILDf response surface...')
      if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
         X = repmat(ILDf_freqaxis',size(dist_ild_signed_sort,2),1);
         Y = ILD_bands_sort';
         [XI,YI] = meshgrid(ILDf_freqaxis,Neuron{cell_num}.tif_ildaxis);
         ILDf_linpredict = griddata(X,Y,ActILD_result2',XI,YI);
         NaN_index = find(isnan(ILDf_linpredict));
         ILDf_ccindex = find(~isnan(ILDf_linpredict));
         ILDf_linpredict(NaN_index) = 0;
      end
   else
      disp('Calculating ILDf response surface...')
      result = sum(result_pieces,1);
      %Calculate modelled ILDf meansurf
      if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
         X = repmat(ILDf_freqaxis',size(dist_ild_signed_sort,2),1);
         Y = ILD_bands_sort';
         [XI,YI] = meshgrid(ILDf_freqaxis,Neuron{cell_num}.tif_ildaxis);
         ILDf_linpredict = griddata(X,Y,ActILD_result',XI,YI);
         NaN_index = find(isnan(ILDf_linpredict));
         ILDf_ccindex = find(~isnan(ILDf_linpredict));
         ILDf_linpredict(NaN_index) = 0;
      end
   end
         
   
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
   
   
   %Find the correlation coefficient between the predicted and measured ILDAlone RS
   temp = corrcoef(IA_linpredict,Neuron{cell_num}.ia_meansurf{1});
   cc(1,cell_num) = temp(1,2)^2;
   
   %Find the correlation coefficient between the Tonal ILD/freq predicted and measured ILDAlone RS
   for cc_loc = 1:length(Neuron{cell_num}.ia_locs{1})
      cc_index(cc_loc) = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,1)...
         & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,2)));
   end
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tia_meanarray')));
      temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},Neuron{cell_num}.tia_meanarray(cc_index));
      cc(2,cell_num) = temp(1,2)^2;
   end
   
  
   %Plot
   mainfig = figure('Units','inches','Position',[0 0 fig_wd fig_ht],'PaperPositionMode','auto');
   
   %plot Measured Tonal ILD/Freq RS
   subplot('Position',[.07 .85 .38 .10]);
   set(gca,'FontSize',8,'FontName','Garamond');
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
      hold on
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
      colormap(cmap);
   end
   plot(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),...
      ILD_matrix_focus(ind_minfreq:ind_maxfreq,IA_bestloc_ind(cell_num)),...
      'yellow','LineWidth',1.5);
   axis([min(ILDf_freqaxis) max(ILDf_freqaxis)...
         min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
   xlabel('Freq (Hz)','FontSize',8);
   ylabel('ILD (dB)','FontSize',8);
   title(['Cell # ' num2str(cell_num) ', Tonal ILD/Freq'],'FontSize',8);
   
   %plot Measured BP ILD/Freq RS
   subplot('Position',[.55 .85 .38 .10]);
   set(gca,'FontSize',8,'FontName','Garamond');
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')))
      plotsurf(Neuron{cell_num}.bif_freqaxis,Neuron{cell_num}.bif_ildaxis,Neuron{cell_num}.bif_meansurf');
      hold on
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
      colormap(cmap);
   end
   plot(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),...
      ILD_matrix_focus(ind_minfreq:ind_maxfreq,IA_bestloc_ind(cell_num)),...
      'yellow','LineWidth',1.5);
   axis([min(ILDf_freqaxis) max(ILDf_freqaxis)...
         min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
   xlabel('Freq (Hz)','FontSize',8);
   ylabel('ILD (dB)','FontSize',8);
   title('BP ILD/Freq','FontSize',8);
   
   %plot Modeled ILD/Freq RS
   subplot('Position',[.30 .65 .40 .10]);
   set(gca,'FontSize',8,'FontName','Garamond');
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
      plotsurf(ILDf_freqaxis,Neuron{cell_num}.tif_ildaxis,ILDf_linpredict);
      hold on
      h = colorbar; set(h,'FontSize',8);
      colormap(cmap);
   end
   plot(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),...
      ILD_matrix_focus(ind_minfreq:ind_maxfreq,IA_bestloc_ind(cell_num)),...
      'yellow','LineWidth',1.5);
   axis([min(ILDf_freqaxis) max(ILDf_freqaxis)...
         min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
   xlabel('Freq (Hz)','FontSize',8);
   ylabel('ILD (dB)','FontSize',8);
   t = Neuron{cell_num}.tif_meansurf';
   temp = corrcoef(t(ILDf_ccindex),ILDf_linpredict(ILDf_ccindex));
   cctonal_ildf(cell_num) = temp(1,2)^2;
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')))
      b = Neuron{cell_num}.bif_meansurf';
      temp = corrcoef(b(ILDf_ccindex),ILDf_linpredict(ILDf_ccindex));
      ccbp_ildf(cell_num) = temp(1,2)^2;
      title(str2mat([Model ' ILD/Freq'],['Tonal r^{2} = ' num2str(cctonal_ildf(cell_num))...
            ', BP r^{2} = ' num2str(ccbp_ildf(cell_num))]),'FontSize',8);
   else
      title(str2mat([Model ' ILD/Freq'],['Tonal r^{2} = ' num2str(cctonal_ildf(cell_num))]),...
         'FontSize',8);
   end
   
   

	%plot True Space RS
   subplot('Position',[.07 .35 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
   	plotdiam(Neuron{cell_num}.ts_azi{1},Neuron{cell_num}.ts_ele{1},Neuron{cell_num}.ts_diamond{1});
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
   	colormap(cmap);
   end
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   title(['True Space'],'FontSize',8);
   if(exist('space_contour'))
      hold on;
      plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
      hold off;
   end
      
   
   %plot ILDAlone RS
   subplot('Position',[.40 .35 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},Neuron{cell_num}.ia_diamond{1});
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
   colormap(cmap);
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   title('Measured ILDAlone','FontSize',8);
   if(exist('space_contour'))
      hold on;
      plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
      hold off;
   end
   
   %plot ILDAlone RS Predicted from Tonal Data
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   subplot('Position',[.70 .35 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   plotdiam(Neuron{cell_num}.tia_azi,Neuron{cell_num}.tia_ele,Neuron{cell_num}.tia_diamond);
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
   colormap(cmap);
   title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc(2,cell_num))]),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
	end

   %plot ILDAlone RS Predicted from Model
   subplot('Position',[.07 .07 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},IA_linpredict_diam);
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
   colormap(cmap);
   title(str2mat([Model ' Predicted ILDAlone'],['r^{2} = ', num2str(cc(1,cell_num))]),'FontSize',8);
   if(exist('space_contour'))
      hold on;
      plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
      hold off;
   end
   
   %plot Measured ILDAlone vs. Tonal ILDAlone
   subplot('Position',[.40 .07 .25 .25])
   set(gca,'FontSize',8,'FontName','Garamond')
   plot(Neuron{cell_num}.ia_meansurf{1}/max(Neuron{cell_num}.ia_meansurf{1}),...
      Neuron{cell_num}.tia_meanarray(cc_index)/max(Neuron{cell_num}.tia_meanarray(cc_index)),'k.')
   hold on
   line(0:.1:1,0:.1:1,'Color','red');
   xlabel('Measured ILDAlone'); ylabel('Line Integral ILDAlone');
   title('Scatter Plot of Line Integral Model');
   axis square
   
   %plot Measured ILDAlone vs. Linmodel9_LR ILDAlone
   subplot('Position',[.73 .07 .25 .25])
   set(gca,'FontSize',8,'FontName','Garamond')
   plot(Neuron{cell_num}.ia_meansurf{1}/max(Neuron{cell_num}.ia_meansurf{1}),...
      IA_linpredict/max(IA_linpredict),'k.')
   hold on
   line(0:.1:1,0:.1:1,'Color','red');
   xlabel('Measured ILDAlone'); ylabel('Linmodel9_{LR} ILDAlone');
   title('Scatter Plot of Linmodel9_{LR} Model');
   axis square
   
   %Print
   if(print_flag)
      figure(firstfig)
      print -loose
      
      figure(secondfig)
      print -loose
      
      figure(mainfig)
      print -loose
   end
   
   end %end if to ask for tif_meansurf
   clear space_contour
end
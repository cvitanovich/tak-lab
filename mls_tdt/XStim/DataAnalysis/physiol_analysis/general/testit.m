%Linmodel9 - script descibing and testing a linear model of spectral integration in ILDAlone RS
clear;
Model = 'Linmodel9_{LR}';

%set parameters
bird_number = 899;
begin_cells = 1;
end_cells = 24;
pca_weights = 0;
temporal_flag = 0;
time_limits = [110 120];
tonal_weights = 0;

%Graphics parameters
cmap = jet;
linecolor = 'white';
bkgndcolor = [0.85 0.85 0.85];
fig_ht = 10;
fig_wd = 8;
Dot_markersize = 4;
plot_dots = 0; %flag for plotting Activity/ILD and Activity/ILD_dist curves
plot_corr = 0; %flag for plotting the correlations

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
clear ITD_matrix ABI_matrix test_numbers reps

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
      %1. ILD_frequency frequency axis, 2. ILDAlone meansurf and diamond
      
      if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
         ILDf_freqaxis = Neuron{cell_num}.bif_freqaxis;
      else
         ILDf_freqaxis = Neuron{cell_num}.tif_freqaxis;
      end
      
      if(temporal_flag == 1)
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
   for num_loc = 1:size(IA_locs,1)
      ILDmat_index(num_loc) =...
         max(find(HRTFinfo.location_matrix(1,:) == IA_locs(num_loc,1) &...
         HRTFinfo.location_matrix(2,:) == IA_locs(num_loc,2)));
   end
   ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   
   %Get the frequency indices
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      ILDf_freqaxis(length(ILDf_freqaxis))));
   end
   
   %Find the best ILDAlone RS location
   if(temporal_flag == 1)
      IA_bestloc_ind = [...
            129   144    14   198   163   214   163   200   282    99   127   250 ...
            130   179   342   159   197    78   218   109     0    95   125   128];
   else
      IA_bestloc_ind(cell_num) = max(find(IA_meansurf == max(IA_meansurf)));
   end
   if(cell_num == 1) IA_bestloc_ind(cell_num) = 129; end
   if(cell_num == 4) IA_bestloc_ind(cell_num) = 198; end
   if(cell_num == 7) IA_bestloc_ind(cell_num) = 163; end
   if(cell_num == 9) IA_bestloc_ind(cell_num) = 282; end
   if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
   if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
   if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
   if(cell_num == 18) IA_bestloc_ind(cell_num) = 78; end
   if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
   
	%TIF data
   [TIA_meanarray2,...
         dirs,full_rs, new_f_ax, new_ild_ax] = if2space4(ILD_matrix',...
      HRTFinfo.hrtf_freqs, HRTFinfo.location_matrix,...
      Neuron{cell_num}.tif_meansurf, Neuron{cell_num}.tif_freqaxis, Neuron{cell_num}.tif_ildaxis);
   [TIA_azi2,TIA_ele2,TIA_diamond2] = ...
      array2diamond(TIA_meanarray2,HRTFinfo.location_matrix);
   
   %Plot
   mainfig = figure('Units','inches','Position',[0 0 fig_wd fig_ht],'PaperPositionMode','auto');
   %plot ILDAlone RS
   subplot('Position',[.07 .35 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   plotdiam(IA_azi,IA_ele,IA_diamond);
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
   subplot('Position',[.40 .35 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   plotdiam(Neuron{cell_num}.tia_azi,Neuron{cell_num}.tia_ele,Neuron{cell_num}.tia_diamond);
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
      colormap(cmap);
   temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},Neuron{cell_num}.tia_meanarray(ILDmat_index));
   cc(2,cell_num) = temp(1,2)^2;   
   title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc(2,cell_num))]),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
	end
   
   %plot ILDAlone RS Predicted from Tonal Data using if2space4
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   subplot('Position',[.70 .35 .25 .25])
   set(gca,'FontSize',8,'Color',bkgndcolor,'FontName','Garamond')
   axis square
   plotdiam(TIA_azi2,TIA_ele2,TIA_diamond2);
      h = colorbar; set(h,'FontSize',8,'FontName','Garamond');
      colormap(cmap);
      
   temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},TIA_meanarray2(ILDmat_index));
   cc(1,cell_num) = temp(1,2)^2;   
   title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc(1,cell_num))]),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
	end
   
end
end
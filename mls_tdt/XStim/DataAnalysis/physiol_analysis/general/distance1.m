%distance1: a script to analyze how a cell's location-dependent ILDAlone activity is a function
%of the distance of a location's ILD spectrum from that of the cell's optimal ILD spectrum,
%for different frequency bands

clear;

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
load 'd:\mlspezio\matlab\save\Neuron'

%get the HRTF info
%[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = ...
%   get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
%clear ITD_matrix ABI_matrix

%Begin looping over the neurons
count = 0;
for cell_num = 6:6
   count = count + 1;
   
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
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - Neuron{cell_num}.tif_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      Neuron{cell_num}.tif_freqaxis(length(Neuron{cell_num}.tif_freqaxis))));
   
   %Find the best ILDAlone RS location
   IA_bestloc_ind = max(find(Neuron{cell_num}.ia_meansurf{1} == max(Neuron{cell_num}.ia_meansurf{1})));
   
   %Define Frequency Bandlimits
   freq_incr = 500;
   count = 0;
   for freq = 2000:500:10500
      count = count + 1;
      freq_bands(count,:) = [freq freq+freq_incr];
   end
   if(isempty(freq_bands))
      size_freq_bands = size(new_if_freqaxis,1);
   else
      size_freq_bands = size(freq_bands,1);
   end
   
   %Compute distance for each location for each frequency band
   for band_num = 1:size(freq_bands,1)
      [y,minfreqindex] = min(abs(HRTFinfo.hrtf_freqs - freq_bands(band_num,1)));
      [y,maxfreqindex] = min(abs(HRTFinfo.hrtf_freqs - freq_bands(band_num,2)));
      ILD_opt = repmat(ILD_matrix_focus(minfreqindex:maxfreqindex,IA_bestloc_ind),...
         1,size(ILD_matrix_focus,2));
      dist(band_num,:) = (sqrt(mean((ILD_matrix_focus(minfreqindex:maxfreqindex,:) -...
         ILD_opt).^2,1)).^0.6);
   end
   
   %Find how many figures will be necessary
   num_standard_figs = 6;
   num_tot_figs = num_standard_figs + size_freq_bands;
   num_colsinplot = 4;
   num_rowsinplot = ceil(num_tot_figs/num_colsinplot);
   
   %Find the correlation coefficient between the Tonal ILD/freq predicted and measured ILDAlone RS
   for cc_loc = 1:length(Neuron{cell_num}.ia_locs{1})
      cc_index(cc_loc) = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,1)...
         & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,2)));
   end
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tia_meanarray')));
   temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},Neuron{cell_num}.tia_meanarray(cc_index));
   cc(2,count) = temp(1,2)^2;
   end
  
   %Plot
   fig_ht = 9.5;
   fig_wd = 7.5;
   mainfig = figure;
   set(mainfig,'Units','inches');
   set(mainfig,'Position',[0 0 fig_wd fig_ht]);
   
   %plot True Space RS
   ts_h = subplot(num_rowsinplot,num_colsinplot,1);
   axis square
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
   plotdiam(Neuron{cell_num}.ts_azi{1},Neuron{cell_num}.ts_ele{1},Neuron{cell_num}.ts_diamond{1});
   colorbar
   colormap(cmap);
   end
   set(ts_h,'Color',bkgndcolor);
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   title(['Cell ' num2str(cell_num) ', True Space'],'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
      
   %plot ILD/Freq RS
   if_h = subplot(num_rowsinplot,num_colsinplot,2);
   hold on
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
   colorbar
   colormap(cmap);
   end
   plot(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),ILD_matrix_focus(ind_minfreq:ind_maxfreq,IA_bestloc_ind),...
      'yellow','LineWidth',1.5);
   axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
         min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
   xlabel('Freq (Hz)','FontSize',8);
   ylabel('ILD (dB)','FontSize',8);
   title('Tonal ILD/Freq','FontSize',8);
   
   %plot ILDAlone RS
   ia_h = subplot(num_rowsinplot,num_colsinplot,4);
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},Neuron{cell_num}.ia_diamond{1});
   colorbar
   colormap(cmap);
   set(ia_h,'Color',bkgndcolor);
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   title('Measured ILDAlone','FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
   norm_IA_surface = Neuron{cell_num}.ia_meansurf{1}./max(Neuron{cell_num}.ia_meansurf{1});
   
   %plot ILDAlone RS Predicted from Tonal Data
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
   tia_h = subplot(num_rowsinplot,num_colsinplot,5);
   axis square   
   plotdiam(Neuron{cell_num}.tia_azi,Neuron{cell_num}.tia_ele,Neuron{cell_num}.tia_diamond);
   colorbar
   colormap(cmap);
   set(tia_h,'Color',bkgndcolor);
   title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc(2,count))]),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
end
end
   %plot the Activity vs. Distance plots for each frequency band
   for band_num = 1:size(freq_bands,1)
      subplot(num_rowsinplot,num_colsinplot,num_standard_figs + band_num);
      plot(dist(band_num,:),norm_IA_surface,'.');
      temp = corrcoef(dist(band_num,:),norm_IA_surface);
      cc_ad(cell_num,band_num) = temp(1,2);
      title(str2mat(['Freq = ' num2str(freq_bands(band_num,1)) '-' num2str(freq_bands(band_num,2))],...
            ['r^{2} = ', num2str(cc_ad(cell_num,band_num))]),'FontSize',8);
   end
      

   clear space_contour
end
%Linmodel9 - script descibing and testing a linear model of spectral integration in ILDAlone RS
clear;
Model = 'Linmodel9';

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

%get the HRTF info
%[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = ...
%   get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
%clear ITD_matrix ABI_matrix

%Begin looping over the neurons
count = 0;
for cell_num = 19:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
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
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - Neuron{cell_num}.tif_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      Neuron{cell_num}.tif_freqaxis(length(Neuron{cell_num}.tif_freqaxis))));
   end
   
   %Find the best ILDAlone RS location
   IA_bestloc_ind = max(find(Neuron{cell_num}.ia_meansurf{1} == max(Neuron{cell_num}.ia_meansurf{1})));
   if(cell_num == 1) IA_bestloc_ind = 129; end
   if(cell_num == 4) IA_bestloc_ind = 198; end
   if(cell_num == 7) IA_bestloc_ind = 163; end
   if(cell_num == 9) IA_bestloc_ind = 282; end
   if(cell_num == 11) IA_bestloc_ind = 127; end
   if(cell_num == 13) IA_bestloc_ind = 130; end
   if(cell_num == 16) IA_bestloc_ind = 159; end
   if(cell_num == 18) IA_bestloc_ind = 78; end
   if(cell_num == 24) IA_bestloc_ind = 128; end
   
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
   
   %minimize the model error
   Options = optimset('Display','iter','MaxIter',8000,'TolFun',0.002,'TolX',0.01,'MaxFunEvals',24000);
   V = [10*ones(1,size(freq_bands,1)) zeros(1,size(freq_bands,1))...
         1]; %ILD tuning widths & ILD peak shifts for model
   [V,fval,exitflag,output] = fminsearch('errfun_specint1',...
      V,...
      Options,...
      Neuron{cell_num}.ia_meansurf{1},...
      new_if_surface,...
      new_if_freqaxis,...
      Neuron{cell_num}.tif_ildaxis,...
      ILD_matrix_focus(ILDmat_freqind(1):ILDmat_freqind(2),:),...
      IA_bestloc_ind,...
      freq_bands);
   
   %Calculate Model prediction with bestfit values V
   
   [err,result,V,freq_weight] = errfun_specint1(V,...
      Neuron{cell_num}.ia_meansurf{1},...
      new_if_surface,...
      new_if_freqaxis,...
      Neuron{cell_num}.tif_ildaxis,...
      ILD_matrix_focus(ILDmat_freqind(1):ILDmat_freqind(2),:),...
      IA_bestloc_ind,...
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
   
   f_weights(cell_num,:) = freq_weight;
   
   %Find the correlation coefficient between the predicted and measured ILDAlone RS
   temp = corrcoef(IA_linpredict,Neuron{cell_num}.ia_meansurf{1});
   cc(1,count) = temp(1,2)^2;
   vars(count,:) = V;
   
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
   ts_h = subplot(3,3,1);
   set(ts_h,'Units','normalized');
   set(ts_h,'Position',[0.1 0.7 0.25 0.25]);
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
   if_h = subplot(3,3,2);
   hold on
   set(if_h,'Units','normalized');
   set(if_h,'Position',[0.45 0.7 0.50 0.25]);
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
   ia_h = subplot(3,3,4);
   set(ia_h,'Units','normalized');
   set(ia_h,'Position',[0.1 0.4 0.25 0.25]);
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
   
   %plot ILDAlone RS Predicted from Tonal Data
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
   tia_h = subplot(3,3,5);
   set(tia_h,'Units','normalized');
   set(tia_h,'Position',[0.4 0.4 0.25 0.25]);
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

   %plot ILDAlone RS Predicted from Model
   mia_h = subplot(3,3,6);
   set(mia_h,'Units','normalized');
   set(mia_h,'Position',[0.7 0.4 0.25 0.25]);
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},IA_linpredict_diam);
   colorbar
   colormap(cmap);
   set(mia_h,'Color',bkgndcolor);
   title(str2mat([Model ' Predicted ILDAlone'],['r^{2} = ', num2str(cc(1,count))]),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
   
   if(length(Neuron{cell_num}.ia_meansurf) > 1)
   %plot ILDAlone RS of frequency band 1
   ia_h = subplot(3,3,7);
   set(ia_h,'Units','normalized');
   set(ia_h,'Position',[0.1 0.1 0.25 0.25]);
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{2},Neuron{cell_num}.ia_ele{2},Neuron{cell_num}.ia_diamond{2});
   colorbar
   colormap(cmap);
   set(ia_h,'Color',bkgndcolor);
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   bandlimits = [Neuron{cell_num}.ia_header{2}(14) - (Neuron{cell_num}.ia_header{2}(15)/2) ...
         Neuron{cell_num}.ia_header{2}(14) + (Neuron{cell_num}.ia_header{2}(15)/2)];
   title(str2mat('Measured ILDAlone',['freq band = (' num2str(bandlimits(1)), ',',...
         num2str(bandlimits(2)), ')' ]),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
   end
   
   if(length(Neuron{cell_num}.ia_meansurf) > 1)
   %plot ILDAlone RS of frequency band 2
   ia_h = subplot(3,3,9);
   set(ia_h,'Units','normalized');
   set(ia_h,'Position',[0.7 0.1 0.25 0.25]);
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{4},Neuron{cell_num}.ia_ele{4},Neuron{cell_num}.ia_diamond{4});
   colorbar
   colormap(cmap);
   set(ia_h,'Color',bkgndcolor);
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   bandlimits = [Neuron{cell_num}.ia_header{4}(14) - (Neuron{cell_num}.ia_header{4}(15)/2) ...
         Neuron{cell_num}.ia_header{4}(14) + (Neuron{cell_num}.ia_header{4}(15)/2)];
   title(str2mat('Measured ILDAlone',['freq band = (' num2str(bandlimits(1)), ',',...
         num2str(bandlimits(2)), ')']),'FontSize',8);
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
   end
   
   end %end if to ask for tif_meansurf
   clear space_contour
end
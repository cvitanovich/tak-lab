%Linmodel1 - script descibing and testing a linear model of spectral integration in ILDAlone RS
clear; %close all

%set parameters
bird_number = 899;
cmap = jet;
bkgndcolor = [0.85 0.85 0.85];
%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;

%specify model
[err,func] = errorfun_spectralintegration3;

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron'

%get the HRTF info
[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = ...
   get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);

%Begin looping over the neurons
count = 0;
for cell_num = 19:20
   count = count + 1;
   
   %Get the contour for the SRF of the neuron
   [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
      srfstats(Neuron{cell_num}.ts_azi{1},...
      Neuron{cell_num}.ts_ele{1},...
      Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
   
   %Fit the model and get the predicted ILDAlone RS
   IA_fullspec = Neuron{cell_num}.ia_meansurf{1};
   for band = 1:2
      IA_bandwidths(band,:) = Neuron{cell_num}.ia_meansurf{2*band};
   end
      IA_locations = Neuron{cell_num}.ia_locs{1};
      %parse the bandwidths into LEFT and RIGHT space
      [index1,IA_piece1] = parse_srs(IA_bandwidths(1,:),IA_locations',[-90 -5],[]);
      [index2,IA_piece2] = parse_srs(IA_bandwidths(1,:),IA_locations',[0 90],[]);
      [index3,IA_piece3] = parse_srs(IA_bandwidths(2,:),IA_locations',[-90 -5],[]);
      [index4,IA_piece4] = parse_srs(IA_bandwidths(2,:),IA_locations',[0 90],[]);
         
   %minimize the model error
   Options = optimset('Display','iter','MaxIter',500,'TolFun',0.002,'MaxFunEvals',1000);
   V = [1 1 1 1]; %coefficients for model
   [V,fval,exitflag,output] = fminsearch('errorfun_spectralintegration5',V,Options,...
      IA_fullspec,IA_piece1,IA_piece2,IA_piece3,IA_piece4,...
      index1,index2,index3,index4);
   
   lofreq = zeros(size(IA_fullspec,1),size(IA_fullspec,2));
   hifreq = zeros(size(IA_fullspec,1),size(IA_fullspec,2));
   lofreq(index1) = V(1) * IA_piece1;
   lofreq(index2) = V(2) * IA_piece2;
   hifreq(index3) = V(3) * IA_piece3;
   hifreq(index4) = V(4) * IA_piece4;
   result = lofreq + hifreq;
   
   IA_linpredict = result;
   
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
   temp = corrcoef(IA_linpredict,IA_fullspec);
   cc(1,count) = temp(1,2)^2;
   vars(count,:) = V;
   
   %Find the correlation coefficient between the Tonal ILD/freq predicted and measured ILDAlone RS
   for cc_loc = 1:length(Neuron{cell_num}.ia_locs{1})
      cc_index(cc_loc) = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,1)...
         & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(cc_loc,2)));
   end
   temp = corrcoef(IA_fullspec,Neuron{cell_num}.tia_meanarray(cc_index));
   cc(2,count) = temp(1,2)^2;
   
   %Find the best ILDAlone RS location
   ind_maxact = max(find(Neuron{cell_num}.ia_meansurf{1} == max(Neuron{cell_num}.ia_meansurf{1})));
   ind_bestloc = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(ind_maxact,1)...
      & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(ind_maxact,2)));
   
   %Get the HRTF frequency indices
   min_freq = min(Neuron{cell_num}.tif_freqaxis);
   max_freq = max(Neuron{cell_num}.tif_freqaxis);
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   
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
   plotdiam(Neuron{cell_num}.ts_azi{1},Neuron{cell_num}.ts_ele{1},Neuron{cell_num}.ts_diamond{1});
   colorbar
   colormap(cmap);
   set(ts_h,'Color',bkgndcolor);
   xlabel('Azimuth (deg)','FontSize',8);
   ylabel('Elevation (deg)','FontSize',8);
   title('True Space, measured','FontSize',8);
   hold on;
   plot(space_contour(1,:),space_contour(2,:),'black','LineWidth',1.5);
   hold off;
   
   %plot ILD/Freq RS
   if_h = subplot(3,3,2);
   hold on
   set(if_h,'Units','normalized');
   set(if_h,'Position',[0.45 0.7 0.50 0.25]);
   plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
   colorbar
   colormap(cmap);
   plot(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),ILD_matrix(ind_minfreq:ind_maxfreq,ind_bestloc),...
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
   hold on;
   plot(space_contour(1,:),space_contour(2,:),'black','LineWidth',1.5);
   hold off;
   
   %plot ILDAlone RS Predicted from Tonal Data
   tia_h = subplot(3,3,5);
   set(tia_h,'Units','normalized');
   set(tia_h,'Position',[0.4 0.4 0.25 0.25]);
   axis square
   plotdiam(Neuron{cell_num}.tia_azi,Neuron{cell_num}.tia_ele,Neuron{cell_num}.tia_diamond);
   colorbar
   colormap(cmap);
   set(tia_h,'Color',bkgndcolor);
   title(str2mat('Tonal Predicted ILDAlone',['r^{2} = ', num2str(cc(2,count))]),'FontSize',8);
   hold on;
   plot(space_contour(1,:),space_contour(2,:),'black','LineWidth',1.5);
   hold off;
   
   %plot ILDAlone RS Predicted from Model
   mia_h = subplot(3,3,6);
   set(mia_h,'Units','normalized');
   set(mia_h,'Position',[0.7 0.4 0.25 0.25]);
   axis square
   plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},IA_linpredict_diam);
   colorbar
   colormap(cmap);
   set(mia_h,'Color',bkgndcolor);
   title(str2mat('Model Predicted ILDAlone',['r^{2} = ', num2str(cc(1,count)) ', coeffs = '],...
      ['(' num2str(vars(count,1)), ' ',...
         num2str(vars(count,2)), ' ',...
         num2str(vars(count,3)), ' ',...
         num2str(vars(count,4)), ')']),'FontSize',8);
   hold on;
   plot(space_contour(1,:),space_contour(2,:),'black','LineWidth',1.5);
   hold off;
   
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
   hold on;
   plot(space_contour(1,:),space_contour(2,:),'black','LineWidth',1.5);
   hold off;
   
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
   hold on;
   plot(space_contour(1,:),space_contour(2,:),'black','LineWidth',1.5);
   hold off;
end

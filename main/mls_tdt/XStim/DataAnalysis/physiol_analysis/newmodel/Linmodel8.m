%Linmodel4 - script descibing and testing a linear model of spectral integration in ILDAlone RS
clear;
Model = 'Linmodel6';

%set parameters
bird_number = 899;
cmap = 1-gray;
linecolor = 'white';
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
clear ITD_matrix ABI_matrix

%Begin looping over the neurons
count = 0;
for cell_num = 19:20
   count = count + 1;
   
   %Get the contour for the SRF of the neuron
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
   [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
      srfstats(Neuron{cell_num}.ts_azi{1},...
      Neuron{cell_num}.ts_ele{1},...
      Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
   end
   
   %Find the best ILDAlone RS location
   ind_maxact = max(find(Neuron{cell_num}.ia_meansurf{1} == max(Neuron{cell_num}.ia_meansurf{1})));
   IA_bestloc = ind_maxact;
   ind_bestloc = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(ind_maxact,1)...
      & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(ind_maxact,2)));
   
   %Get the HRTF frequency indices
   min_freq = min(Neuron{cell_num}.tif_freqaxis);
   max_freq = max(Neuron{cell_num}.tif_freqaxis);
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece1
   min_freq = 2000; max_freq = 3000;
   [y,ind_minfreq_1] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_1] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece2
   min_freq = 3000; max_freq = 4000;
   [y,ind_minfreq_2] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_2] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece3
   min_freq = 4000; max_freq = 5000;
   [y,ind_minfreq_3] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_3] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece4
   min_freq = 5000; max_freq = 6000;
   [y,ind_minfreq_4] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_4] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece5
   min_freq = 6000; max_freq = 7000;
   [y,ind_minfreq_5] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_5] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece6
   min_freq = 7000; max_freq = 8000;
   [y,ind_minfreq_6] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_6] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece7
   min_freq = 8000; max_freq = 9000;
   [y,ind_minfreq_7] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_7] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   %ILD_piece8
   min_freq = 9000; max_freq = 10000;
   [y,ind_minfreq_8] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq_8] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
     
   %Get the HRTF location indices
   IA_locs = Neuron{cell_num}.ia_locs{1}';
   for num_loc = 1:size(IA_locs,2)
      loc_index(num_loc) = ...
         max(find(HRTFinfo.location_matrix(1,:) == IA_locs(1,num_loc) & ...
         HRTFinfo.location_matrix(2,:) == IA_locs(2,num_loc)));
   end
   
   ILD_piece1 = ILD_matrix(ind_minfreq_1:ind_maxfreq_1,loc_index);
   ILD_piece2 = ILD_matrix(ind_minfreq_2:ind_maxfreq_2,loc_index);
   ILD_piece3 = ILD_matrix(ind_minfreq_3:ind_maxfreq_3,loc_index);
   ILD_piece4 = ILD_matrix(ind_minfreq_4:ind_maxfreq_4,loc_index);
   ILD_piece5 = ILD_matrix(ind_minfreq_5:ind_maxfreq_5,loc_index);
   ILD_piece6 = ILD_matrix(ind_minfreq_6:ind_maxfreq_6,loc_index);
   ILD_piece7 = ILD_matrix(ind_minfreq_7:ind_maxfreq_7,loc_index);
   ILD_piece8 = ILD_matrix(ind_minfreq_8:ind_maxfreq_8,loc_index);
   
   %Fit the model and get the predicted ILDAlone RS
   IA_fullspec = Neuron{cell_num}.ia_meansurf{1};
   
   %minimize the model error
   Options = optimset('Display','iter','MaxIter',2000,'TolFun',0.002,'TolX',0.01,'MaxFunEvals',4000);
   V = [ones(1,8) ...
         10*ones(1,8)...
         zeros(1,size(ILD_piece1,1))...
         zeros(1,size(ILD_piece2,1))...
         zeros(1,size(ILD_piece3,1))...
         zeros(1,size(ILD_piece4,1))...
         zeros(1,size(ILD_piece5,1))...
         zeros(1,size(ILD_piece6,1))...
         zeros(1,size(ILD_piece7,1))...
         zeros(1,size(ILD_piece8,1))...
      ]; %Frequency weights, ILD tuning widths, ILD peak shifts
   [V,fval,exitflag,output] = fminsearch('errorfun_spectralintegration14',...
      V,...
      Options,...
      IA_fullspec,...
      ILD_piece1,...
      ILD_piece2,...
      ILD_piece3,...
      ILD_piece4,...
      ILD_piece5,...
      ILD_piece6,...
      ILD_piece7,...
      ILD_piece8,...
      IA_bestloc);
   
   %Calculate Model prediction with bestfit values V
max_activity = max(IA_fullspec)/8;

%specify indexing for V
start_weights = 1; end_weights = 8;
start_ILDwidth = 9; end_ILDwidth = 16;
start_ILDpeakshift1 = 17; end_ILDpeakshift1 = start_ILDpeakshift1 + size(ILD_piece1,1)-1;
start_ILDpeakshift2 = end_ILDpeakshift1 + 1; end_ILDpeakshift2 = start_ILDpeakshift2 + size(ILD_piece2,1)-1;
start_ILDpeakshift3 = end_ILDpeakshift2 + 1; end_ILDpeakshift3 = start_ILDpeakshift3 + size(ILD_piece3,1)-1;
start_ILDpeakshift4 = end_ILDpeakshift3 + 1; end_ILDpeakshift4 = start_ILDpeakshift4 + size(ILD_piece4,1)-1;
start_ILDpeakshift5 = end_ILDpeakshift4 + 1; end_ILDpeakshift5 = start_ILDpeakshift5 + size(ILD_piece5,1)-1;
start_ILDpeakshift6 = end_ILDpeakshift5 + 1; end_ILDpeakshift6 = start_ILDpeakshift6 + size(ILD_piece6,1)-1;
start_ILDpeakshift7 = end_ILDpeakshift6 + 1; end_ILDpeakshift7 = start_ILDpeakshift7 + size(ILD_piece7,1)-1;
start_ILDpeakshift8 = end_ILDpeakshift7 + 1; end_ILDpeakshift8 = start_ILDpeakshift8 + size(ILD_piece8,1)-1;

%IA_piece1
V_mat = repmat(V(start_ILDpeakshift1:end_ILDpeakshift1)',1,size(ILD_piece1,2));
ILD_opt_piece1 = repmat(ILD_piece1(:,IA_bestloc),1,size(ILD_piece1,2));
diff_ILD_piece1 = sqrt(mean((ILD_piece1 - ILD_opt_piece1 + V_mat).^2,1));
result_piece1 = max_activity * (1 - (diff_ILD_piece1/V(start_ILDwidth)));

%IA_piece2
V_mat = repmat(V(start_ILDpeakshift2:end_ILDpeakshift2)',1,size(ILD_piece2,2));
ILD_opt_piece2 = repmat(ILD_piece2(:,IA_bestloc),1,size(ILD_piece2,2));
diff_ILD_piece2 = sqrt(mean((ILD_piece2 - ILD_opt_piece2 + V_mat).^2,1));
result_piece2 = max_activity * (1 - (diff_ILD_piece2/V(start_ILDwidth+1)));

%IA_piece3
V_mat = repmat(V(start_ILDpeakshift3:end_ILDpeakshift3)',1,size(ILD_piece2,2));
ILD_opt_piece3 = repmat(ILD_piece3(:,IA_bestloc),1,size(ILD_piece3,2));
diff_ILD_piece3 = sqrt(mean((ILD_piece3 - ILD_opt_piece3 + V_mat).^2,1));
result_piece3 = max_activity * (1 - (diff_ILD_piece3/V(start_ILDwidth+2)));

%IA_piece4
V_mat = repmat(V(start_ILDpeakshift4:end_ILDpeakshift4)',1,size(ILD_piece2,2));
ILD_opt_piece4 = repmat(ILD_piece4(:,IA_bestloc),1,size(ILD_piece4,2));
diff_ILD_piece4 = sqrt(mean((ILD_piece4 - ILD_opt_piece4 + V(12)).^2,1));
result_piece4 = max_activity * (1 - (diff_ILD_piece4/V(start_ILDwidth+3)));

%IA_piece5
V_mat = repmat(V(start_ILDpeakshift5:end_ILDpeakshift5)',1,size(ILD_piece2,2));
ILD_opt_piece5 = repmat(ILD_piece5(:,IA_bestloc),1,size(ILD_piece5,2));
diff_ILD_piece5 = sqrt(mean((ILD_piece5 - ILD_opt_piece5 + V_mat).^2,1));
result_piece5 = max_activity * (1 - (diff_ILD_piece5/V(start_ILDwidth+4)));

%IA_piece6
V_mat = repmat(V(start_ILDpeakshift6:end_ILDpeakshift6)',1,size(ILD_piece2,2));
ILD_opt_piece6 = repmat(ILD_piece6(:,IA_bestloc),1,size(ILD_piece6,2));
diff_ILD_piece6 = sqrt(mean((ILD_piece6 - ILD_opt_piece6 + V_mat).^2,1));
result_piece6 = max_activity * (1 - (diff_ILD_piece6/V(start_ILDwidth+5)));

%IA_piece7
V_mat = repmat(V(start_ILDpeakshift7:end_ILDpeakshift7)',1,size(ILD_piece2,2));
ILD_opt_piece7 = repmat(ILD_piece7(:,IA_bestloc),1,size(ILD_piece7,2));
diff_ILD_piece7 = sqrt(mean((ILD_piece7 - ILD_opt_piece7 + V(15)).^2,1));
result_piece7 = max_activity * (1 - (diff_ILD_piece7/V(start_ILDwidth+6)));

%IA_piece8
V_mat = repmat(V(start_ILDpeakshift8:end_ILDpeakshift8)',1,size(ILD_piece2,2));
ILD_opt_piece8 = repmat(ILD_piece8(:,IA_bestloc),1,size(ILD_piece8,2));
diff_ILD_piece8 = sqrt(mean((ILD_piece8 - ILD_opt_piece8 + V_mat).^2,1));
result_piece8 = max_activity * (1 - (diff_ILD_piece8/V(start_ILDwidth+7)));


result = (V(1)*result_piece1) + (V(2)*result_piece2)...
   + (V(3)*result_piece3) + (V(4)*result_piece4) + (V(5)*result_piece5)...
   + (V(6)*result_piece6) + (V(7)*result_piece7) + (V(8)*result_piece8);

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
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
   end
   
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
   if(exist('space_contour'))
   hold on;
   plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
   hold off;
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
   title(str2mat([Model ' Predicted ILDAlone'],['r^{2} = ', num2str(cc(1,count))],['ILD width = (', ...
         num2str(vars(count,17)), ')']),'FontSize',8);
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
   
   clear space_contour
end

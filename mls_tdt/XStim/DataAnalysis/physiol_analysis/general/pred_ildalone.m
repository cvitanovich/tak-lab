function [Tonal_testpars,...
   	Tonal_azi,...
      Tonal_ele,...
      Tonal_diamond,...
      BP_azi,...
      BP_ele,...
      BP_diamond,...
      ILDAlone_azi,...
      ILDAlone_ele,...
      ILDAlone_diamond,...
      ITD_meancurve,...
      ITD_diamond,...
      ABI_curve,...
      ABI_diamond,...
      Tonal_cc,...
      BP_cc...
   ] = pred_ildalone(bird_number,side_of_brain,test_numbers,get_HRTF,get_ITDmatrix);
%This script produces the following plots for a particular recording site
%	1. ILD alone data
%	2. ILD alone predicted from tonal and 1/3-BP noise ILD/Freq surfaces
%	3. Tonal ILD/Freq surface
%	4. 1/3-BP noise ILD/Freq surface

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];

%Get information from bird's HRTF catalogue - for use in predicting ILD Alone Surfaces
%Header information
[filetype,...
info_blocks,...
n_channels,n_lines, sampling_rate,...
first_line,last_line,num_locations,...
comment1, comment2] = mtlrh([hrtf_directory hrtf_file]);
%Location information
temp = mtlrdir([hrtf_directory hrtf_file]);
location_matrix = sph2dbl(temp);
itddir = location_matrix;
clear temp

hrtf_freqs = [first_line:(last_line-first_line)/(n_lines-1):last_line]';

if (get_HRTF == 1)
	%Transfer Functions
	for loc = 1:num_locations
   	temp = mtlrch([hrtf_directory hrtf_file],2*loc-1);
      left_raw(:,loc) = temp;
      left = 20*log10(abs(temp));
      clear temp
      
   	temp = mtlrch([hrtf_directory hrtf_file],2*loc);
      right_raw(:,loc) = temp;
      right = 20*log10(abs(temp));
      clear temp


		ILD_matrix(:,loc) = right - left;
		ABI_matrix(:,loc) = right + left ./2;
      
      if(mod(loc,10) == 0)
      	disp(['Finished location ' num2str(loc)])
   	end
   end
   eval(['save ' hrtf_directory num2str(bird_number) 'ildmatrix ILD_matrix left_raw right_raw']);
   eval(['save ' hrtf_directory num2str(bird_number) 'abimatrix ABI_matrix']);
else
   eval(['load ' hrtf_directory num2str(bird_number) 'ildmatrix.mat']);
   eval(['load ' hrtf_directory num2str(bird_number) 'abimatrix.mat']);
end

if(get_ITDmatrix == 1)
   [IPD_matrix,hrtf_ipd_freqs] = make_ipdmatrix(hrtf_directory,hrtf_file2,2000,11500);
   [ITD_matrix,hrtf_itd_freqs,azi2,ele2] = sunwrap_matrix(IPD_matrix,hrtf_ipd_freqs,location_matrix);
	eval(['save ' hrtf_directory 'ipdmatrix IPD_matrix hrtf_ipd_freqs']);
   eval(['save ' hrtf_directory 'itdmatrix ITD_matrix hrtf_itd_freqs']);
else
 	eval(['load ' hrtf_directory 'ipdmatrix.mat']);
   eval(['load ' hrtf_directory 'itdmatrix.mat']);
end

%Get the mean response surface of the Tonal ILD/Freq data
[Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(1), 1, 0);

		%Input to proc_test:    proc_test899(bird, side, test, cluster, pltflg, count_starttime, ...
      %								count_starttime, count_endtime, include_reps);
      
%Get the mean response surface of the BP ILD/Freq data
[BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(2), 1, 0);

%Get the mean response surface of the ILD Alone data
[ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(3), 1, 0);
[ILDAlone_azi,ILDAlone_ele,temp] = array2diamond(ILDAlone_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(ILDAlone_azi, ILDAlone_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(locs(:,1)==el & locs(:,2)==az))
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
	
   ILDAlone_diamond = intval + temp;
   ILDAlone_diamond = ILDAlone_diamond./max(max(ILDAlone_diamond));
   
if (test_numbers(4) ~= -1) %check to see if we actually have a True Space measurement
%Get the mean response surface of the True Space data
[Space_meansurf, Space_stdsurf, Space_dim2vals, Space_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(4), 1, 0);
[Space_azi,Space_ele,temp] = array2diamond(Space_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Space_azi, Space_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(locs(:,1)==el & locs(:,2)==az))
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
	
   Space_diamond = intval + temp;
   Space_diamond = Space_diamond./max(max(Space_diamond));
 end  
 
%Generate the ITD Spatial diamond 
%1. Get the ITD curve (measured using Broadband Noise)
[ITD_meancurve, ITD_stdcurve, ITD_dim2vals, itdfreqs, ITD_testpars, ...
      spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(5), 1, 0);
%2. Get the Mean Frequency Response Curve from the Tonal ILD/Freq surface
Tonal_freqcurve = max(Tonal_meansurf,[],2);
%3. Get the ITD diamond 
[ITD_diamond, ITD_az, ITD_el, wa_itd] = itd2space3(ITD_meancurve,...
   itdfreqs, ITD_matrix', itddir, hrtf_itd_freqs', Tonal_freqcurve, Tonal_dim1vals);
ITD_diamond = ITD_diamond./max(max(ITD_diamond));

%Generate the ABI Spatial diamond
%1. Get the ABI curve (simulated using "sigmoid" function)
ABI_axis = -200:0;
ABI_curve = sigmoid(ABI_axis,1.3,test_numbers(6)); %values betw -1 and 0
%2. Get the ABI diamond using the abi2space function
[ABI_diamond,ABI_el,ABI_az] = abi2space(ABI_matrix,hrtf_freqs,location_matrix,...
   ABI_axis,ABI_curve,Tonal_freqcurve,Tonal_dim1vals,test_numbers(6));
   
%Generate the PREDICTED ILD Alone surface from the Tonal Data
[Tonal_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
   hrtf_freqs, location_matrix, Tonal_meansurf, Tonal_dim1vals, Tonal_dim2vals);
[Tonal_azi,Tonal_ele,Tonal_diamond] = array2diamond(Tonal_prediction,location_matrix);
Tonal_diamond = Tonal_diamond./max(max(Tonal_diamond));

%Generate the PREDICTED True Space surface from the Tonal Data
Tonal_truespace_diamond = Tonal_diamond .* ITD_diamond .* ABI_diamond.^2;
Tonal_truespace_diamond = Tonal_truespace_diamond./max(max(Tonal_truespace_diamond));

%Generate the PREDICTED ILD Alone surface from the BP Data
[BP_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
   hrtf_freqs, location_matrix, BP_meansurf, BP_dim1vals, BP_dim2vals);
[BP_azi,BP_ele,BP_diamond] = array2diamond(BP_prediction,location_matrix);
BP_diamond = BP_diamond./max(max(BP_diamond));

%Generate the PREDICTED True Space surface from the BP Data
BP_truespace_diamond = BP_diamond .* ITD_diamond .* ABI_diamond;
BP_truespace_diamond = BP_truespace_diamond./max(max(BP_truespace_diamond));

%Plot the relevant activity surfaces
set(0,'DefaultAxesCreateFcn','set(gca,''Color'',''blue'')');
allplot_h = figure;
set(allplot_h,'units','normal');
set(allplot_h,'Position',[0 0 0.8 0.8]);

%Tonal ILD/Freq surface
tonalsurf_h = subplot(3,3,1);
set(tonalsurf_h,'Position',[0.05 0.86 0.3 0.1]);
plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf');
colorbar
title([num2str(bird_number) side_of_brain ' Depth= ' num2str(Tonal_testpars(7)) ' ITD= '...
      num2str(Tonal_testpars(2)) ' ABI= ' num2str(Tonal_testpars(4))])
xlabel('Frequency (Hz)');
ylabel('ILD (dB)');

%BP ILD/Freq surface
BPsurf_h = subplot(3,3,3);
set(BPsurf_h,'Position',[0.55 0.86 0.3 0.1]);
plotsurf(BP_dim1vals,BP_dim2vals,BP_meansurf');
colorbar
xlabel('Frequency (Hz)');

%Tonal ILDAlone diamond
tonaldiam_h = subplot(3,3,4);
set(tonaldiam_h,'Position',[0.05 0.45 0.3 0.3]);
plotdiam(Tonal_azi,Tonal_ele,Tonal_diamond);
set(gca,'Color','none');
title('Pure Tones')

%BP ILDAlone diamond
BPdiam_h = subplot(3,3,5);
set(BPdiam_h,'Position',[0.37 0.45 0.3 0.3]);
plotdiam(BP_azi,BP_ele,BP_diamond);
set(gca,'Color','none');
title('BP Noise')

%Measured ILDAlone diamond
ILDAlonediam_h = subplot(3,3,6);
set(ILDAlonediam_h,'Position',[0.7 0.45 0.3 0.3]);
plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diamond);
set(gca,'Color','none');
title('ILD Alone SRF')

%Tonal Truespace diamond
Tonal_truespace_h = subplot(3,3,7);
set(Tonal_truespace_h,'Position',[0.05 0.05 0.3 0.3]);
plotdiam(Tonal_azi,Tonal_ele,Tonal_truespace_diamond);
set(gca,'Color','none');
title('Tonal Truespace SRF')

%BP Truespace diamond
BP_truespace_h = subplot(3,3,8);
set(BP_truespace_h,'Position',[0.37 0.05 0.3 0.3]);
plotdiam(BP_azi,BP_ele,BP_truespace_diamond);
set(gca,'Color','none');
title('BP Truespace SRF')

%Measured Truespace diamond
if(test_numbers(4) ~= -1)
   Spacediam_h = subplot(3,3,9);
   set(Spacediam_h,'Position',[0.7 0.05 0.3 0.3]);
   plotdiam(Space_azi,Space_ele,Space_diamond);
   set(gca,'Color','none');
   title('True Space RF');
end

%Statistical analysis of comparison

%1. ILDAlone -- Correlation across all measured locations
%Use only those locations that were measured
index = missmask == 0;

%Only use that portion of Tonal & BP prediction that matches the ILDAlone (i.e., w/o ele=90)
Tonal_diamond_2 = Tonal_diamond(1:size(ILDAlone_diamond,1),:);
BP_diamond_2 = BP_diamond(1:size(ILDAlone_diamond,1),:);

cc = corrcoef(ILDAlone_diamond(index),Tonal_diamond_2(index));
Tonal_cc(1,1) = cc(1,2);
clear cc

cc = corrcoef(ILDAlone_diamond(index),BP_diamond_2(index));
BP_cc(1,1) = cc(1,2);
clear cc

%2. Truespace -- Correlation across all measured locations
if(test_numbers(4) ~= -1)
   Tonal_truespace_2 = Tonal_truespace_diamond(1:size(Space_diamond,1),:);
   BP_truespace_2 = BP_truespace_diamond(1:size(Space_diamond,1),:);
   cc = corrcoef(Space_diamond(index),Tonal_truespace_2(index));
   Tonal_cc(1,2) = cc(1,2);
   clear cc
   
   cc = corrcoef(Space_diamond(index),BP_truespace_2(index));
   BP_cc(1,2) = cc(1,2);
   clear cc
end

return


%Script to analyze the impact on Spatial Response Surface (SRS)
%of limiting the frequency bandwidth of a sound
clear;close all

%Tests for analysis
%Each row is a different cell
%Numbers in the following order: [0-13kHz, 3-7kHz, 8-11kHz]
test_numbers = [
   207 212 213
   215 221 222];

%Set parameters
bird_number = 899;
side_of_brain = 'r';
get_hrtf = 0;

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

for cell_num = 1:size(test_numbers,1)

%1. Get the mean response surface of the True Space data with 0-13 kHz
[Space0013_meansurf, Space0013_stdsurf, Space0013_dim2vals, Space0013_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,1), 1, 0);
[Space0013_azi,Space0013_ele,temp] = array2diamond(Space0013_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Space0013_azi, Space0013_ele);

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
	
   Space0013_diamond = intval + temp;
   max_value = max(max(Space0013_diamond));
   Space0013_diamond = Space0013_diamond./max(max(Space0013_diamond));
 
 %2. Get the mean response surface of the True Space data with 3-7 kHz
[Space0307_meansurf, Space0307_stdsurf, Space0307_dim2vals, Space0307_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,2), 1, 0);
[Space0307_azi,Space0307_ele,temp] = array2diamond(Space0307_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Space0307_azi, Space0307_ele);

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
	
   Space0307_diamond = intval + temp;
   Space0307_diamond = Space0307_diamond./max_value;
 
 %3. Get the mean response surface of the True Space data with 8-11 kHz
[Space0811_meansurf, Space0811_stdsurf, Space0811_dim2vals, Space0811_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,3), 1, 0);
[Space0811_azi,Space0811_ele,temp] = array2diamond(Space0811_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Space0811_azi, Space0811_ele);

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
	
   Space0811_diamond = intval + temp;
   Space0811_diamond = Space0811_diamond./max_value;

temp = corrcoef(Space0307_meansurf,Space0013_meansurf);
cc_lof = temp(1,2); clear temp;
temp = corrcoef(Space0811_meansurf,Space0013_meansurf);
cc_hif = temp(1,2); clear temp;

 figure
 subplot(3,1,1)
 plotdiam(Space0307_azi,Space0307_ele,Space0307_diamond);
 set(gca,'Color','none');
 title(['SRS for 3-7 kHz Bandpassed Noise, r = ' num2str(cc_lof)])
 colorbar
 
 subplot(3,1,2)
 plotdiam(Space0811_azi,Space0811_ele,Space0811_diamond);
 set(gca,'Color','none');
 title(['SRS for 8-11 kHz Bandpassed Noise, r = ' num2str(cc_hif)])
 colorbar
 
 subplot(3,1,3)
 plotdiam(Space0013_azi,Space0013_ele,Space0013_diamond);
 set(gca,'Color','none');
 title('SRS for 0-13 kHz Bandpassed Noise')
 colorbar

end %end loop over cells
 

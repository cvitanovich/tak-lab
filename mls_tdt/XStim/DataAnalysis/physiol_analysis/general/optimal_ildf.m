function [Pred_ILDAlone,...
      Pred_ILDAlone_diamond,...
      Pred_ILDAlone_azi,...
      Pred_ILDAlone_ele,...
      optimal_ildf_surface,...
      freq_axis,...
      ild_axis] = optimal_ildf(bird_number,side_of_brain,test_numbers,hrtf_file,get_hrtf)
%Function to generate an ILD/Frequency surface from a cell's *measured* ILDAlone surface
%test_numbers: [Tonal_ILD_Freq ILDAlone_SRF]
%Overall view of function:
%1. Get the bird's HRTF catalogue
%2. Get the cell's measured Tonal ILD/Freq surface & calculate frequency tuning curve from it
%3. Get the cell's measured ILDAlone response surface
%4. Trim the ILD and frequency axes of the HRTF spectra to match those of the ILDAlone RS
%5. Generate the optimal ILD/Freq surface

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];

%1. Get information from bird's HRTF catalogue - for use in predicting ILD Alone Surfaces
   %Header information
   [filetype,...
         info_blocks,...
         n_channels,n_lines, sampling_rate,...
         first_line,last_line,num_locations,...
         comment1, comment2] = mtlrh([hrtf_directory hrtf_file]);
   %Location information
   temp = mtlrdir([hrtf_directory hrtf_file]);
   location_matrix = sph2dbl(temp);
   clear temp
   
   hrtf_freqs = [first_line:(last_line-first_line)/(n_lines-1):last_line]';
   
if(get_hrtf == 1)
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

%2. Get the mean response surface of the Tonal ILD/Freq data
[Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(1), 1, 0);
%freq_tuning = max(Tonal_meansurf,[],2);
%freq_tuning = freq_tuning./max(freq_tuning);
freq_tuning = ones(1,length(Tonal_dim1vals));


%3. Get the mean response surface of the ILD Alone data
[ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(2), 1, 0);
max_ILDAlone_meansurf = max(ILDAlone_meansurf);
   
%4. Trim the ILD and Frequency axes of the HRTF spectra
%maxild = max(Tonal_dim2vals);
%minild = min(Tonal_dim2vals);
maxild = max(max(ILD_matrix));
minild = min(min(ILD_matrix));
maxfreq = max(Tonal_dim1vals);
minfreq = min(Tonal_dim1vals);

maxild = round(maxild*2)/2;
if (maxild>max(ILDAlone_dim2vals)) maxild = maxild-.5; end;  % insures that we didn't round up and out of bounds
minild = round(minild*2)/2;
if (minild<min(ILDAlone_dim2vals)) minild = minild+.5; end;

%trim ILD values
ILD_matrix(ILD_matrix>=maxild) = maxild*ones(size(ILD_matrix(ILD_matrix>=maxild))); %set all values outside range equal to max or min
ILD_matrix(ILD_matrix<=minild) = minild*ones(size(ILD_matrix(ILD_matrix<=minild)));
ILD_matrix = round(ILD_matrix*2)/2;  % limit to .5 db steps
new_hrtf_ild = minild:.5:maxild;
%ild_axis = resample(new_hrtf_ild,length(Tonal_dim2vals),length(new_hrtf_ild));
ild_axis = Tonal_dim2vals;

% trim frequency information
infreqs = hrtf_freqs<=maxfreq & hrtf_freqs>=minfreq;
ILD_matrix = ILD_matrix(infreqs,:);
new_hrtf_freq = hrtf_freqs(infreqs);
freq_axis = Tonal_dim1vals;

%5. Generate the optimal ILDf surface
optimal_ildf_surface_mat = zeros(length(freq_axis),length(ild_axis),length(locs));
nonzero_index = ones(length(freq_axis),length(ild_axis));
count = 0;
for l = 1:length(locs)
%for l = 1:1
   temp = zeros(length(freq_axis),length(ild_axis));
   locindex = max(find(location_matrix(1,:) == locs(l,1) & location_matrix(2,:) == locs(l,2)));
   for findex = 1:length(freq_axis)
      freqindex = max(find(new_hrtf_freq <= freq_axis(findex)));
      if(isempty(freqindex))
         freqindex = min(find(new_hrtf_freq >= freq_axis(findex)));
      end
      ildindex = max(find(ild_axis <= ILD_matrix(freqindex,locindex)));
      if(isempty(ildindex))
         ildindex = min(find(ild_axis >= ILD_matrix(freqindex,locindex)));
      end
      count = count+1;
      ildindex_vec(count) = ildindex;
      temp(findex,ildindex) = 1 * freq_tuning(findex);      
   end
   temp = temp .* ILDAlone_meansurf(l);
   nonzero_index = nonzero_index + (temp ~= 0);
   optimal_ildf_surface_mat(:,:,l) = optimal_ildf_surface_mat(:,:,l) + temp;
   clear temp
end

optimal_ildf_surface = mean(optimal_ildf_surface_mat,3) ./ (nonzero_index);
%optimal_ildf_surface = optimal_ildf_surface - mean(optimal_ildf_surface(:));
optimal_ildf_surface = optimal_ildf_surface./max(max(optimal_ildf_surface));

%Try to generate a predicted ILDAlone space response from the optimal ILD/Freq surface
for l = 1:length(locs)
   temp = zeros(length(freq_axis),length(ild_axis));
   locindex = max(find(location_matrix(1,:) == locs(l,1) & location_matrix(2,:) == locs(l,2)));
   for findex = 1:length(freq_axis)
      freqindex = max(find(new_hrtf_freq <= freq_axis(findex)));
      if(isempty(freqindex))
         freqindex = min(find(new_hrtf_freq >= freq_axis(findex)));
      end
      ildindex = max(find(ild_axis <= ILD_matrix(freqindex,locindex)));
      if(isempty(ildindex))
         ildindex = min(find(ild_axis >= ILD_matrix(freqindex,locindex)));
      end
      temp(findex,ildindex) = 1;
   end
   Pred_ILDAlone(l) = sum(max(optimal_ildf_surface_mat(:,:,l) .* temp,[],2));
   clear temp
end

[Pred_ILDAlone_azi,Pred_ILDAlone_ele,temp] = array2diamond(Pred_ILDAlone,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(Pred_ILDAlone_azi, Pred_ILDAlone_ele);

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
	
   Pred_ILDAlone_diamond = intval + temp;
   Pred_ILDAlone_diamond = Pred_ILDAlone_diamond./max(max(Pred_ILDAlone_diamond));


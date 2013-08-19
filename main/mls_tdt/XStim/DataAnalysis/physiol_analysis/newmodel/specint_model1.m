%A script to model the spectral integration of ICCls and ICx neurons,
%based on the linear model by Nelken, Kim and Young (J Neurophysiol 78(1997)800-811)
%Uses the matlab function fminsearch, which in turn uses the error function
%errorfun_spectralintegration1 (sum squared errors in all ILD/Freq bins)

clear; close all
%Input needed for functions below
bird_number = 899;
side_of_brain = 'r';
test_numbers = [
   76 77 79 80 75 -62
   92 93 94 95 91 -72
   98 100 101 -1 97 -80
   105 106 107 -1 104 -72
   110 111 112 113 109 -75
   123 124 125 126 122 -75
   129 130 131 132 127 -65
   135 136 137 138 134 -73
   143 144 145 146 142 -73
   162 158 159 161 153 -75
   167 168 169 170 166 -70
   173 174 175 176 172 -75];
colormap_var = 'jet';
plotflag = 0;

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF        = 0;

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


for cell_num = 1:size(test_numbers,1)
%1. Get the mean response surface of the Tonal ILD/Freq data
[Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,1), 1, 0);

%2. Get the mean response surface of the BP ILD/Freq data
[BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,2), 1, 0);

%3. Get the frequency range
min_freq_ind = min(find(hrtf_freqs >= min(Tonal_dim1vals)));
max_freq_ind = max(find(hrtf_freqs <= max(Tonal_dim1vals)));
new_freqs = hrtf_freqs(min_freq_ind:max_freq_ind);

%3. Interpolate in frequency
for ILD_num = 1:length(Tonal_dim2vals)
   Tonal_meansurf_2(:,ILD_num) = interp1(Tonal_dim1vals,Tonal_meansurf(:,ILD_num),new_freqs);
   
   BP_meansurf_2(:,ILD_num) = interp1(BP_dim1vals,BP_meansurf(:,ILD_num),new_freqs);
end
disp('Finished interpolating in frequency.')

%4. Minimize the appropriate function
Options = optimset('Display','iter','MaxIter',10,'TolFun',0.02);
var1_start = 200; %effective bandwidth in Hz
var2_start = 3; %exponent for Tonal response; assume linearity first
%V = [var1_start var2_start];
V = [var2_start];
[vars(cell_num),fval,exitflag,output] = fminsearch('errorfun_spectralintegration1',V,Options,...
   new_freqs,Tonal_meansurf_2,BP_meansurf_2);
disp(['Solution for cell ' num2str(cell_num) ' = ' num2str(vars)])

%5. Calculate resultant modeled surface
target_bandwidth = 1/3;
%determine weights
f_weights(1) = 1; %from paper
f_weights(length(new_freqs)) = 1;
for k = 2:(length(new_freqs) - 1)
   f_weights(k) = (new_freqs(k+1) - new_freqs(k-1)) / 2;
end
f_weights = repmat(f_weights',1,size(Tonal_meansurf,2));


for num_freq = 1:length(new_freqs)
   [lo_limit,hi_limit] = bandlimits(new_freqs(num_freq),target_bandwidth);
   [temp,index_1] = min(abs(new_freqs - lo_limit));
   [temp,index_2] = min(abs(new_freqs - hi_limit));
   
   %calculate the new ILD/Freq surface from the Tonal data
   Tonal_meansurf_new(num_freq,:) = ...
      mean( f_weights(index_1:index_2,:).^vars(cell_num) .* (Tonal_meansurf_2(index_1:index_2,:)),1);
end
for num_freq = 1:length(Tonal_dim1vals)
   [temp,index_1] = min(abs(new_freqs - Tonal_dim1vals(num_freq)));
   Tonal_meansurf_resamp(num_freq,:) = Tonal_meansurf_new(index_1,:);
   BP_meansurf_resamp(num_freq,:) = BP_meansurf_2(index_1,:);
end

%6. Display the resultant meansurfaces

figure
subplot(3,1,1)
plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf')
colorbar
title('Tonal Measured')

subplot(3,1,2)
plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf_resamp')
colorbar
title('BP Modeled')

subplot(3,1,3)
plotsurf(Tonal_dim1vals,Tonal_dim2vals,BP_meansurf_resamp')
colorbar
title('BP Measured')

clear f_weights
end %end loop over cells
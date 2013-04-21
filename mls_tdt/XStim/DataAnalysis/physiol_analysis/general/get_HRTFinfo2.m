function [L_matrix,R_matrix,ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = get_HRTFinfo2(bird_number,...
   hrtf_file,...
   hrtf_file2,...
   get_HRTF,...
   get_ITDmatrix)

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];

%Get information from bird's HRTF catalogue - for use in predicting ILD Alone Surfaces
%Header information
[HRTFinfo.filetype,...
HRTFinfo.headerinfo_blocks,...
HRTFinfo.n_channels,HRTFinfo.n_lines, HRTFinfo.sampling_rate,...
HRTFinfo.first_line,HRTFinfo.last_line,HRTFinfo.num_locations,...
HRTFinfo.comment1, HRTFinfo.comment2] = mtlrh([hrtf_directory hrtf_file]);
%Location information
temp = mtlrdir([hrtf_directory hrtf_file]);
HRTFinfo.location_matrix = sph2dbl(temp);
clear temp

HRTFinfo.hrtf_freqs = [HRTFinfo.first_line:(HRTFinfo.last_line-HRTFinfo.first_line)/...
      (HRTFinfo.n_lines-1):HRTFinfo.last_line]';

if (get_HRTF == 1)
	%Transfer Functions
	for loc = 1:HRTFinfo.num_locations
   	temp = mtlrch([hrtf_directory hrtf_file],2*loc-1);
      left_raw(:,loc) = temp;
      left = 20*log10(abs(temp));
      clear temp
      
   	temp = mtlrch([hrtf_directory hrtf_file],2*loc);
      right_raw(:,loc) = temp;
      right = 20*log10(abs(temp));
      clear temp
      
      L_matrix(:,loc) = left;
		R_matrix(:,loc) = right;
		ILD_matrix(:,loc) = right - left;
		ABI_matrix(:,loc) = right + left ./2;
      
      if(mod(loc,10) == 0)
      	disp(['Finished location ' num2str(loc)])
   	end
   end
   eval(['save ' hrtf_directory num2str(bird_number) 'ildmatrix ILD_matrix left_raw right_raw']);
   eval(['save ' hrtf_directory num2str(bird_number) 'abimatrix ABI_matrix']);
   eval(['save ' hrtf_directory num2str(bird_number) 'LR L_matrix R_matrix']);
else
   eval(['load ' hrtf_directory num2str(bird_number) 'ildmatrix.mat']);
   eval(['load ' hrtf_directory num2str(bird_number) 'abimatrix.mat']);
   eval(['load ' hrtf_directory num2str(bird_number) 'LR']);
end

if(get_ITDmatrix == 1)
   [IPD_matrix,hrtf_ipd_freqs] = make_ipdmatrix(hrtf_directory,hrtf_file2,2000,11500);
   [ITD_matrix,hrtf_itd_freqs,azi2,ele2] = sunwrap_matrix(IPD_matrix,hrtf_ipd_freqs,HRTFinfo.location_matrix);
	eval(['save ' hrtf_directory 'ipdmatrix IPD_matrix hrtf_ipd_freqs']);
   eval(['save ' hrtf_directory 'itdmatrix ITD_matrix hrtf_itd_freqs']);
else
 	eval(['load ' hrtf_directory 'ipdmatrix.mat']);
   eval(['load ' hrtf_directory 'itdmatrix.mat']);
end

%Construct ILD_matrix for a given bird's HRTFs

HRTF_file = 'e:\spezio\matlab\scripts\tdt\XStim\HRTFfilts\894f_tdAD';
junkfile = 'junk';
mtlfft(HRTF_file,junkfile,2048);
[filetype,info_blocks,...
      n_channels,n_lines, sampling_rate,...
      first_line,last_line,n_directions,...
      comment1, comment2] = mtlrh(junkfile);

nyq = sampling_rate/2;
hrtf_freqs = 0:nyq/(n_lines-1):nyq;

for loc = 1:n_directions
   L = 20*log10(abs(mtlrch(junkfile,2*loc-1)));
   R = 20*log10(abs(mtlrch(junkfile,2*loc)));
   
   ILD_matrix(:,loc) = R-L;
   disp(['Finished with location ' num2str(loc)])
end

HRTFinfo.hrtf_freqs = hrtf_freqs;
HRTFinfo.location_matrix = mtlrdir(HRTF_file);
save 'e:\spezio\matlab\scripts\tdt\XStim\HRTFfilts\894_ILD_matrix' ILD_matrix HRTFinfo

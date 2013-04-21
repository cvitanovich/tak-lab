%Testhyp1
%Script to test an hypothesis about the relative strengths of contribution of certain
%ILD/freq ranges
clear;

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

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'

%get the HRTF info
[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = ...
   get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
clear ITD_matrix ABI_matrix

%Begin looping over the neurons
count = 0;
for cell_num = 22:22
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
   ind_bestloc = max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(ind_maxact,1)...
      & HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(ind_maxact,2)));
   
   %Get the HRTF frequency indices
   min_freq = min(Neuron{cell_num}.tif_freqaxis);
   max_freq = max(Neuron{cell_num}.tif_freqaxis);
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - min_freq));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - max_freq));
   
   %Get the ILD spectra in lo frequency range
   freqlimits = [2000 6500];
   if(cell_num == 22)
      AZlimits = [-60 -20];
      ELlimits = [-70 -35];
   elseif(cell_num == 20)
      AZlimits = [50 60];
      ELlimits = [-40 -15];
   end
   [freqindex_lo,locindex,lo_ILDmat] = ...
      parse_HRTF(ILD_matrix,HRTFinfo.hrtf_freqs,HRTFinfo.location_matrix,freqlimits,AZlimits,ELlimits);
   
   %Get the ILD spectra in hi frequency range
   freqlimits = [6500 11000];
   [freqindex_hi,locindex,hi_ILDmat] = ...
      parse_HRTF(ILD_matrix,HRTFinfo.hrtf_freqs,HRTFinfo.location_matrix,freqlimits,AZlimits,ELlimits);
   
   figure
   %plot ILD/Freq RS
   hold on
   plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
   colorbar
   colormap(cmap);
   axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
         min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
   xlabel('Freq (Hz)','FontSize',8);
   ylabel('ILD (dB)','FontSize',8);
   title(str2mat('Tonal ILD/Freq',['AZlim = ('...
         num2str(AZlimits(1)) ',' num2str(AZlimits(2)) ')'],...
      ['ELlim = (' num2str(ELlimits(1)) ',' num2str(ELlimits(2)) ')']), 'FontSize',8);
   
   %plot lo freq ILDmat
   plot(HRTFinfo.hrtf_freqs(freqindex_lo),lo_ILDmat,'cyan');
   
   %plot hi freq ILDmat
   plot(HRTFinfo.hrtf_freqs(freqindex_hi),hi_ILDmat,'magenta');
   
   %plot best TrueSpace location
   plot(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),ILD_matrix(ind_minfreq:ind_maxfreq,ind_bestloc),...
      'yellow','LineWidth',1.5);

end %end loop over cells


   


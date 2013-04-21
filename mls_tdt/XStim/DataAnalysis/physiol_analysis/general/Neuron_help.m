%The following are the properties of the structure Neuron
%  ml_position    - medial/lateral pos
%  ap_position    - anterior/posterior ps
%  depth
%  threshhold(noise)
%  sr_testpars    - Spontaneous rate test parameters
%  sr_mean        - Spontaneuos rate mean (scalar)
%  sr_std         - Spontaneous rate standard deviation (scalar)
%  abi_testpars   - Avg Binaural Intensity test parameters
%  abi_abiaxis    - ABI x-axis
%  abi_mean       - ABI mean curve
%  abi_std        - ABI standard deviation
%  itd_testpars   - Interaural Time Difference test parameters
%  itd_itdaxis    - ITD x-axis
%  itd_mean       - ITD mean curve
%  itd_std        - ITD standard deviation
%  abif_testpars  - ABI/Frequency test parameters
%  abif_freqaxis  - ABI/Frequency freq axis
%  abif_abiaxis   - ABI/Frequency ABI axis
%  abif_meansurf  - ABI/F mean response surface
%  abif_stdsurf   - ABI/F standard deviation of response surface
%  tif_testpars   - Tonal ILD/Freq test parameters
%  tif_freqaxis   - TIF frequency axis
%  tif_ildaxis    - TIF ILD axis
%  tif_meansurf   - TIF mean response surface
%  tif_stdsurf    - TIF standard deviation of response surface
%  bif_testpars   - BandPass ILD/Freq test parameters
%  bif_freqaxis   - BIF frequency axis
%  bif_ildaxis    - BIF ILD axis
%  bif_meansurf   - BIF mean response surface
%  bif_stdsurf    - BIF standard deviation of response surface
%ILDAlone - these properties each have multiple entries due to the different frequency ranges
%  ia_testpars    - ILDAlone test parameters
%  ia_locs        - ILDAlone test locations
%  ia_azi         - Azimuth for ILDAlone diamond
%  ia_ele         - Elevation for ILDAlone diamond
%  ia_meanarray   - ILDAlone mean response surface in array format
%  ia_stdarray    - ILDAlone standard deviation of response surface
%  ia_diamond     - ILDAlone diamond RS
%  tia_meanarray  - PREDICTED Tonal ILDAlone array
%  tia_azi        - Azimuth for PREDICTED Tonal ILDAlone diamond RS
%  tia_ele        - Elevation for PREDICTED Tonal ILDAlone diamond RS
%  tia_diamond    - PREDICTED Tonal ILDAlone diamond RS
%  bia_meanarray  - PREDICTED BP ILDAlone array
%  bia_azi        - Azimuth for PREDICTED BP ILDAlone diamond RS
%  bia_ele        - Elevation for PREDICTED BP ILDAlone diamond RS
%  bia_diamond    - PREDICTED BP ILDAlone diamond RS
%True Space - these properties each have multiple entries due to the different frequency ranges
%  ts_testpars    - True Space test parameters
%  ts_locs        - True Space test locations
%  ts_meanarray   - True Space mean response surface in array format
%  ts_stdarray    - True Space standard deviation of response surface
%  ts_diamond     - True Space diamond
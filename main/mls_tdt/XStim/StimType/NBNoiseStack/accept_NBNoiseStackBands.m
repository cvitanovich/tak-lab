%accept_NBNoiseStackBands: script to accept the specified frequency bands for the NBNoise Stack

XStimParams.spec_NBNoiseStackBands = [];

for bandnum = 1:XStimParams.num_NBNoiseStackBands
   XStimParams.spec_NBNoiseStackBands(bandnum,:) = (sscanf(get(H.spec_NBNoiseStackBands(bandnum),'String'),'%i'))';
end

close(H.NBNoiseStackBands)
H.NBNoiseStackBands = [];
H.num_NBNoiseStackBands = [];

return
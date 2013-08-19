function T = make_targs(Neuron)
%Requires tif and ia psth's

temp_tif = Neuron.tif_psth(:);
numlocs = length(Neuron.ia_psth);

for timept = 1:length(Neuron.psth_bins)
   for meas_num = 1:length(Neuron.ia_psth) + length(temp_tif)
      if(meas_num <= length(Neuron.ia_psth))
         T{timept}(meas_num) = Neuron.ia_psth{meas_num}(timept);
      end
      if(meas_num > length(Neuron.ia_psth))
         T{timept}(meas_num) = temp_tif{meas_num - numlocs}(timept);
      end
   end
end
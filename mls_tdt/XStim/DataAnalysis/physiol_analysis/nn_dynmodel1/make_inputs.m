function P = make_inputs(Neuron,HRTFinfo,L_mat,R_mat,bandwidth)
%P is a cell matrix of matrices. Each cell is a time point, each column in the matrix
%is an input vector p. First row = first input; Second row = second input.
%Number of rows in the plain matrix = number of elements in the input vector.
%Number of cols in the plain matrix = number of concurrent inputs (replicants)

numlocs = length(Neuron.ia_psth);

min_Lmat = min(min(L_mat));
min_Rmat = min(min(R_mat));
L_mat = (L_mat + abs(min_Lmat))/abs(min_Lmat);
R_mat = (R_mat + abs(min_Rmat))/abs(min_Rmat);

for loc_num = 1:length(Neuron.ia_locs{1})
   loc_index(loc_num) = min(find(HRTFinfo.location_matrix(1,:) == Neuron.ia_locs{1}(loc_num,1) &...
      HRTFinfo.location_matrix(2,:) == Neuron.ia_locs{1}(loc_num,2)));
end
L_mat_new = L_mat(:,loc_index);
R_mat_new = R_mat(:,loc_index);

[freq_mat,ild_mat] = meshgrid(Neuron.tif_freqaxis,Neuron.tif_ildaxis);
ild_mat = ild_mat'/max(Neuron.tif_ildaxis);
temp_tif = Neuron.tif_psth(:);

threshhold = 20; %dB
if(bandwidth == 0)
   
      for meas_num = 1:length(Neuron.ia_psth) + length(temp_tif)
         freqs = Neuron.tif_freqaxis;
         for freq_num = 1:length(freqs)
            if(meas_num <= length(Neuron.ia_psth))
               [y,freqind] = min(abs(HRTFinfo.hrtf_freqs - freqs(freq_num)));
               P{1,1}(freq_index,meas_num) = L_mat_new(freq_num,meas_num);
               P{2,1}(freq_num,meas_num) = R_mat_new(freq_num,meas_num);
            end
            
            if(meas_num > length(Neuron.ia_psth))
               [freq_index,ild_index] = ind2sub(size(Neuron.tif_meansurf),meas_num - numlocs);
               ildval = ild_mat(freq_num,ild_index);
               if(ildval < 0)
                  P{1,1}(freq_index,meas_num) = abs(ildval);
                  P{2,1}(freq_index,meas_num) = 0;
               end
               if(ildval > 0)
                  P{1,1}(freq_index,meas_num) = 0;
                  P{2,1}(freq_index,meas_num) = abs(ildval);
               end
               if(ildval == 0)
                  P{1,1}(freq_index,meas_num) = 1;
                  P{2,1}(freq_index,meas_num) = 1;
               end
            end
         end
      end
      
      
else
   
      for meas_num = 1:length(Neuron.ia_psth) + length(temp_tif)
         freqs = Neuron.tif_freqaxis;
         for freq_num = 1:length(freqs)
            [minfreq,maxfreq] = bandlimits(freqs(freq_num),bandwidth);
            
            if(meas_num <= length(Neuron.ia_psth))
               [y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
               [y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
               P{1,1}(freq_num,meas_num) = ...
                  mean(L_mat_new(minfreqind:maxfreqind,meas_num));
               P{2,1}(freq_num,meas_num) = ...
                  mean(L_mat_new(minfreqind:maxfreqind,meas_num));
            end
            
            if(meas_num > length(Neuron.ia_psth))
               [freq_index,ild_index] = ind2sub(size(Neuron.tif_meansurf),meas_num - numlocs);
               ildval = ild_mat(freq_num,ild_index);
               if(ildval <= 0)
                  P{1,1}(freq_index,meas_num) = abs(ildval);
                  P{2,1}(freq_index,meas_num) = 0;
               end
               if(ildval > 0)
                  P{1,1}(freq_index,meas_num) = 0;
                  P{2,1}(freq_index,meas_num) = abs(ildval);
               end
               if(ildval == 0)
                  P{1,1}(freq_index,meas_num) = 1;
                  P{2,1}(freq_index,meas_num) = 1;
               end
            end
         end
      end
   
end

P = repmat(P,1,length(Neuron.psth_bins));

return
      
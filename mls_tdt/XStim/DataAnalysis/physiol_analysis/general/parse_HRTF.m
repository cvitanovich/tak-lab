function [freqindex,locindex,parsed_result] = ...
   parse_HRTF(HRTF_matrix,frequencies,locations,freqlimits,AZlimits,ELlimits);
%Parses an HRTF matrix according to LOCATION criteria
%The matrix must have frequencies assigned to rows, locations assigned to columns
%Usage:
%   parse_HRTF(HRTF_matrix,locations,freqlimits,AZlimits,ELlimits);
%freqlimits: frequency range; AZlimits: azimuth range; ELlimits: elevation range
%Can leave either freqlimits, AZlimits or ELlimits empty to include all values in range

   %Select frequency range
   if(~isempty(freqlimits))
      minfreq = freqlimits(1); maxfreq = freqlimits(2);
      [y,findex_min] = min(abs(frequencies - minfreq));
      [y,findex_max] = min(abs(frequencies - maxfreq));
      freqindex = findex_min:findex_max; 
      temp_HRTFmat = HRTF_matrix(freqindex,:);
   else
      temp_HRTFmat = HRTF_matrix;
   end
   
   if(~isempty(AZlimits) & ~isempty(ELlimits))
      minAZ = AZlimits(1); maxAZ = AZlimits(2); minEL = ELlimits(1); maxEL = ELlimits(2);
      index = find(locations(1,:) >= minEL & locations(1,:) <= maxEL ...
         & locations(2,:) >= minAZ & locations(2,:) <= maxAZ);
   elseif(~isempty(AZlimits) & isempty(ELlimits))
      minAZ = AZlimits(1); maxAZ = AZlimits(2);
      index = find(locations(2,:) >= minAZ & locations(2,:) <= maxAZ);
   elseif(isempty(AZlimits) & ~isempty(ELlimits))
      minEL = ELlimits(1); maxEL = ELlimits(2);
      index = find(locations(1,:) >= minEL & locations(1,:) <= maxEL);
   end
   
   locindex = index;
   parsed_result = temp_HRTFmat(:,locindex);
   return


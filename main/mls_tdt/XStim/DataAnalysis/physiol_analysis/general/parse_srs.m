function [index,parsed_result] = ...
   parse_srs(Spatial_Response_Surface,locations,AZlimits,ELlimits);
%Parses a spatial mean response surface (an array) according to LOCATION criteria
%Usage:
%   parse_srs(Spatial_Response_Surface,spatial_locations,[minAZ maxAZ],[minEL maxEL]);
%minAZ : minimum azimuth, etc.
%Can leave either AZ or EL empty to include all values in range

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

parsed_result = Spatial_Response_Surface(index);

return
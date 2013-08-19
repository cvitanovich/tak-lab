function [resp, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ildspechrtf, specfreqs, locs, ildfreq_data, f_ax, ild_ax, checklocs); 
% IF2SPACE       function which transforms an ild-freq response into spatial tuning curve
% [resp, dirs, full_rs, new_f_ax, new_ild_ax] = 
%                if2space(ildspechrtf, specfreqs, locs, ildfreq_data, f_ax, ild_ax, checklocs); 
%
%                checklocs and locs dimension is  [2 x number of locations]
%                ildspechrtf if the hrtf spectral data (locs along rows, freq long cols)
%                specfreqs is the axis for the hrtf spectrum
%                locs is the direction matrix which comes in the hrtf file [el; az] 
%                ildfreq_data is the response surface for the ild-freq test
	
if (nargin<7) 
  checklocs = locs; 
  loci = 1:size(checklocs,2);
else
  for i = 1:size(checklocs,2)
     loci(i) = find(locs(1, :)==checklocs(1,i) & locs(2,:)==checklocs(2,i));
     if isempty(loci(i)) error('mismatch between locations specified and locations in spectrum'); end;
  end;
end;
						
maxild = max(ild_ax);
minild = min(ild_ax);

minfreq = min(f_ax);
maxfreq = max(f_ax);

% limit ild values in hrtf spectrum to range given in the ild-freq plot
% esssentially this insures that ild values outside of range are counted 
% at the last known ild value in the pure tone data

maxild = round(maxild*2)/2;
if (maxild>max(ild_ax)) maxild = maxild-.5; end;  % insures that we didn't round up and out of bounds
minild = round(minild*2)/2;
if (minild<min(ild_ax)) minild = minild+.5; end;

ildspechrtf(ildspechrtf>=maxild) = maxild*ones(size(ildspechrtf(ildspechrtf>=maxild)));
ildspechrtf(ildspechrtf<=minild) = minild*ones(size(ildspechrtf(ildspechrtf<=minild)));
ildspechrtf = round(ildspechrtf*2)/2;  % limit to .5 db steps
ild_axis_hrtf = minild:.5:maxild;

% trim frequencies from ild hrtf spectrum
[y,minfreqind] = min(abs(specfreqs - minfreq));
[y,maxfreqind] = min(abs(specfreqs - maxfreq));
ildspechrtf = ildspechrtf(:,minfreqind:maxfreqind);
freq_axis_hrtf = specfreqs(minfreqind:maxfreqind);

% interpolate missing ild-freq vals
full_rs = interp2(ild_ax, f_ax, ildfreq_data, ild_axis_hrtf, freq_axis_hrtf, 'cubic');
full_rs(isnan(full_rs)) = 0; %added 9/00
new_f_ax = freq_axis_hrtf;
new_ild_ax = ild_axis_hrtf;

% now take each location and compute neural response
ildinums = ildspechrtf/.5+find(ild_axis_hrtf==0);
nfreqs = size(freq_axis_hrtf,1);
freqinums = 1:nfreqs;

resp = zeros(1,size(checklocs,2));
for i=1:size(checklocs,2)     % for each spatial location
	curspec = ildinums(loci(i),:);
	inums = (curspec-1)*nfreqs+freqinums;
	resp(i) = sum(full_rs(inums))/nfreqs;
end;
		
dirs = checklocs;

return;

% the long-hand way to do it, as a check (above algorithm checks ok)
for i=1:size(checklocs,2)     % for each spatial location
	curspec = ildinums(loci(i),:);
	cursum = 0;
	for j = 1:nfreqs
		curildinum = curspec(j);
	    cursum = cursum + full_rs(j,curildinum);
	end;
	resp(i) = cursum/nfreqs;

end;

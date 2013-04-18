function sample=mksample(fg,N,samplerate,fftlen);
% MKSAMPLE  Return nonuniformly spaced sampling points.
%	    Syntax: MKSAMPLE(fg,spfr,samplerate,fftlen)
%	    fg specifies the upper cutoff-frequencies of different
%	    frequency groups. MKSAMPLE tries do determine N integer
%	    sampling points per frequency group.

fr=samplerate/fftlen;		% frequency resolution
if size(fg,1)>1 fg=fg'; end;	% fg must be column vector
fg=[0 fg];			
for i=2:length(fg)
  sfd=(fg(i)-fg(i-1))/N;	% sample frequency distance in current frequency group
  sffg=(0:(N-1))*sfd+fg(i-1);	% sampling frequencies in current frequency group
  sf=[sf sffg];			% all sampling frequencies
end;
sample=zeros(1,length(sf));
j=2;
for i=2:length(sf)		% search integer k such that k*fr is close to sf(i)
  k=sf(i)/fr;
  k=fix(k);
  if k>sample(j-1)  		% ensure that all samples are unique
    sample(j)=k;
    j=j+1;
  end;
end;
sample=sample(1:(j-1))+1;
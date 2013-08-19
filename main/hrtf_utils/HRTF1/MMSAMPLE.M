function sample=mmsample(H,minsample,maxsample);
% MMSAMPLE Compute sampling points near the maxima and minima points
%          of log10(abs(H)).
%          sample=MMSAMPLE(H,minsample,maxsample);
%          sample does not contain any values below MINSAMPLE or
%          above MAXSAMPLE.

Mprime=derive(smooth(log10(abs(H)),1),1);

sample=find(sign([Mprime(:);1])~=sign([1;Mprime(:)]));
sample=sample(1:length(sample)-1);
sample=sample(find(sample<maxsample & sample>minsample));
sample=sample(:)';
if sample(1)>minsample sample=[minsample sample]; end;

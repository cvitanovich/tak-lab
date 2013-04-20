function [sample,samplestr]=getsample(mode)
% GETSAMPLE return sampling points for for monaural
%           and interaural optimization.
%           [sample,samplestr]=getsample(mode)
%           mode=='m' returns samples for the 
%           optimization of the monaural transfer-
%           function while mode=='i' returns samples 
%           for the optimization of the interaural 
%           transfer function.

if nargin<1 mode='m'; end;

if mode(1)=='m'
   samplestr='[10:2:49 50:4:199 200:15:499 500:20:929 930:50:1400]';
elseif mode(1)=='i'
  samplestr='[10:10:49 50:40:199 200:5:499 500:10:929 930:20:1400]';
else
  error('invalid mode');
end;

eval(sprintf('sample=%s;',samplestr));


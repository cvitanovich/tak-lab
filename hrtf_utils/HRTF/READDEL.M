function [d,mindelay]=readdel(catname)
% READDEL  returns the delay times and the minimum delay time of the
%          impulse responses in the specified catalogue.
%          The first row contains the delay times for the left channels,
%          the second row the time for the right channels
%          [d,mindelay]=READDEL(catname)

if exist([catname , '.mat'])~=2
  error('catalogue not found'); 
end;
eval(sprintf('load %s',catname));
d=delay;
mindelay=min(min(d));ÿ

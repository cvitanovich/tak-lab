function [d,ndirs]=readdirs(catname)
% READDIRS returns the available directions and the number of
%          directions stored in the specified catalogue.
%          The first column of the (2*ndirs) matrix d contains
%          the azimuth angles, the second column the elevation angles.
%          [d,ndirs]=READDIRS(catname)

if exist([catname , '.mat'])~=2
  error('catalogue not found');
end;
eval(sprintf('load %s',catname));
d=flt(:,1:2);
ndirs=size(d,1);

function [locs, data_array] = diam2array(a, e, diam, desired_locs);
% DIAM2ARRAY unwraps a diamond into an 1-D array
% 		[locs data_array] = diam2array(a, e, diam, desired_locs);
% a and e are standard 1D axes for elev and azim
% diam is the numerical values in square matrix representing space
% desired_locs is optional and lets you specify which [el az] locations to
% extract data from
% desired locs, like the locs output argument, is a 2 by number-of-locations
% array where first row is elevation and second row is azimuth
% if desired locs are given, output locs is just the same, otherwise data
% is provided for every location 

[A E] = meshgrid(a,e);

if nargin<4
	data_array = diam(~isnan(diam));  	% this will automatically unwrap, provided there are
													% NaN's in the data

	locs(1,:) = E(~isnan(diam))';
	locs(2,:) = A(~isnan(diam))';
else
	for i = 1:size(desired_locs,2)
		data_array(i) = diam(A==desired_locs(2,i) & E==desired_locs(1,i));
		locs = desired_locs;	
	end;
end;

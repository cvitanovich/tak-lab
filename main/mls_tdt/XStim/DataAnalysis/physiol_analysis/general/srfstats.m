function [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...
            srfstats(az, el, spacediam, limval, plotflag); 
% SRFSTATS spatial receptive field stats
%
% [rf_center_az rf_center_el srflim rf_area, rf_i, space_contour] = ...
%            srfstats(az, el, spacediam, limval, plotflag); 
%
%          extracts the following parameters from a space test
%          rf_center_az   azimuth of SRF center
%          rf_center_el   elevation of SRF center
%          srflim. {lo up l r}  structure with upper, lower, left and right limits of SRF
%          rf_area        number of locations covered by SRF
%          rf_i           index to all locations which were counted in SRF
%          space_contour   points needed to plot the contour used
%                         use:
%									 plot(space_contour(1,:),space_contour(2,:), 'w');
%          
%          limval is the cut-off value for inclusion in the SRF
%                 it is the ratio relative to peak (e.g., .5)
%          computation is based on the main peak only.  Side peaks are ignored.
%          center is found via wieghted mean over rf center area
%    have not figured out how to get a contour without plotting it, so this
%    function will draw over existing plot even if plotflag is not set

if (nargin<4) limval = .5; end;
if (nargin<5) plotflag = 0; end;

norm_diam = spacediam/max(max(spacediam));

[space_contour] = contourc(az,el,norm_diam,[limval limval]);

% deal with cases where there are is more than one contour
% find peak and pick the contour which surrounds it
[A E] = meshgrid(az,el);
maxr = max(max(norm_diam));
baz = A(norm_diam==maxr);
bel = E(norm_diam==maxr);
bestloc = [baz(1) bel(1)];

peakok = 0;
start_i = 1;
while (peakok~=1)
   c1 = space_contour(:,start_i+1:start_i+space_contour(2,start_i)); % needed because first pair is not data, its index
           peakok = min(c1(1,:))<=bestloc(1) & max(c1(1,:))>=bestloc(1) & ...
                    min(c1(2,:))<=bestloc(2) & max(c1(2,:))>=bestloc(2);
           start_i = start_i + space_contour(2,start_i)+1;
end;
space_contour = c1;

% now find limits
srflim.lo = min(space_contour(2,:)); 
srflim.up = max(space_contour(2,:));
srflim.l  = min(space_contour(1,:));
srflim.r  = max(space_contour(1,:));

% do weighted sum to find rf center
keepi = A<=(srflim.r) & ...
        A>=(srflim.l) & ...		
        E<=(srflim.up) & ...
        E>=(srflim.lo);
rf_all = norm_diam>limval;
rf_i = rf_all & keepi;
rf_area = sum(sum(rf_i));

center_vals = norm_diam(rf_i);
%center_vals(isnan(center_vals)) = zeros(size(center_vals(isnan(center_vals))));
center_az = A(rf_i);
center_el = E(rf_i);
rf_center_az = sum(center_vals.*center_az)/sum(center_vals);
rf_center_el = sum(center_vals.*center_el)/sum(center_vals);


if (plotflag)
   plotdiam(az,el,spacediam);
   hold on;
   plot(space_contour(1,:),space_contour(2,:), 'w');
   hold off;
end;

contour_pts = space_contour(1:2,:)';

return;

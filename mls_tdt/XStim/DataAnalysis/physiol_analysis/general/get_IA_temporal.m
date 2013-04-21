function [IA_meansurf,IA_stdsurf,IA_locs,IA_diamond,IA_azi,IA_ele] = ...
   get_IA_temporal(bird_number,side_of_brain,test_number,cluster,plotflag,time_limits)
%Function to return ILDAlone data restricted to a temporal range
%bird_number:				code of bird
%side_of_brain:			l or r
%test_number:				number of ILDAlone test
%time_limits:				[begin end]



[IA_meansurf,...
      IA_stdsurf,...
      dim2vals,...
      dim1vals,...
      testpars,spont_spikes, spont_dur, nincl_reps,...
      IA_locs,...
      header] = proc_test899(bird_number,...
      side_of_brain,...
      test_number,...
      cluster,...
      plotflag,...
      time_limits(1),...
      time_limits(2));
   
[IA_azi,IA_ele,temp] = ...
         array2diamond(IA_meansurf,IA_locs');
      %Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(IA_azi,IA_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(IA_locs(:,1)==el & IA_locs(:,2)==az))
	       missmask(AZ==az & EL==el) = 1;
	    else
	       missmask(AZ==az & EL==el) = 0;
	    end; 
	  end;
	end;
	
	% replace all valid locations presently containing NaN's with zeros
	ri = (missmask==1);
	temp(ri) = zeros(size(temp(ri)));
	% generate interpolation function - uses surrounding 4 squares to generate missing square
   intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
	temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
	intval_old = conv2(temp, intfun,'same')/4;
	intval = intval_old.*missmask;
	
   IA_diamond = intval + temp;
   
return

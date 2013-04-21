function [mean_resp_surf,...
      std_surf,...
      dim2vals,...
      dim1vals,...
      testpars,...
      spont_spikes,...
      spont_dur,...
      nincl_reps,...
      locs,...
      test_head] = proc_test899(isi,...
   
% PROC_TEST897	 process a test and plot results
%
%[resp_surf, std_surf, xaxis, yaxis, testpar, spont_spikes, spont_dur, nincl_reps] = proc_test890(bird, ...
%																	side, ...
%																test, ...
%																	cluster, ... 
%																plot_flag, ...
%																	count_starttime, ...
%																	count_endtime, ...
%															include_reps);
%	cluster 
%					array of cluster numbers to include in the analysis (e.g. [1 2])
%	count_starttime, count_endtime (optional)
%					ms times bracketting the spikes to
%				sum for final count.
%					default is 0 to stim_dur
%	include_reps
%					repetition blocks to include in the final count.
%					default is all reps
%   testpars     is   [freq itd ild abi medlat antpost depth stimdur] or NaN if that dimension varied
%
%	spont_spikes = number of spikes during null periods before and after test
%   spont_dur    = length of time over which spont_spikes were sampled
%   nincl_reps     = the number of reps tested (includes half-finished blocks)


% note: we can use set(1, 'paperposition', [3 .5 3 10]) to set the raster display
% to display better [left bottom width height]

% fill in missing values not supplied by user
if (~exist('pltflg', 'var')) 
	plot_flag = 1;
else
	plot_flag = pltflg;
end;

 
% test types: 1=ILD, 2=ITD, 3=PHASE, 4=FREQ, 5=ABI, 6=ILD-ITD
% 7 = ILD-FREQ, 8 = ILD-ABI, 9 = ITD-FREQ, 10 = FREQ-ABI, 11 = PHASE-ABI, 12 = RATT-LATT,
% 13 = ILD-ITD-FREQ

% set some constants

rast_overrun = min(isi, 50);
if (isi==0) rast_overrun = 20; end;
intrinsic_latency = 5;  % this may be due to apparatus, etc.  All spikes are shifted this much
                        % before display and analysis (RC tests are handled separately, though);
                        % I picked this number based on response offset latency
                        % in fact, it cooresponds well with the latency between trigger onset and the 
                        % beginning of stimulus play
int_latency_space = 21; % same thing, intrinsic latency, but for space.  The delay between trigger
                        % and stimulus onset is more like 23 ms for space stimuli
default_countstart = 0; % seems to be better to always start counting a little after stim begins


% why are there two offsets?  intrinsic latency is due to the stimulus delay.  It takes that long 
% for the stimulus to start during a trial.  countstart is a measure of how long it takes the cell 
% to start firing once it hears the sound.  You can tell the intrinsic latency often by looking at 
% when the stimulus offset is.  


if (~exist('spf_flag'))  % this code will be obsolete when all cache files contain this variable
   spf_flag = 0;
end;
if (~exist('bg_flag'))  % this code will be obsolete when all cache files contain this variable
   bg_flag = 0;
end;

% set start and end of the count period

if (~exist('count_starttime', 'var') | isempty(count_starttime)) 
	if (~rc_flag)
	  start_count = default_countstart;   % this is used for counting
	  count_starttime = default_countstart;  % this is used to display lines on rasters
	else
	  start_count = 0;
	  count_starttime = [];
	end;
	if (bg_flag)
       start_count = abs(bg_lead)+default_countstart;
       count_starttime = start_count;
    end;
    if (use_testinfo & ~isempty(curtestinfo.start))
       start_count = curtestinfo.start;
       count_starttime = curtestinfo.start;
    end;
else
	start_count = count_starttime;
end;
if (~exist('count_endtime', 'var') | isempty(count_endtime))
	end_count = stim_dur + start_count;
	count_endtime = [];
	if (bg_flag)
       end_count = stim_dur;
       count_endtime = end_count;
    end;
    if (use_testinfo & ~isempty(curtestinfo.end))
       end_count = curtestinfo.end;
       count_endtime = curtestinfo.end;
    end;
else
	end_count = count_endtime;
end;	 
disp(['Spikes counted from ' num2str(start_count) ' ms to ' num2str(end_count) ' ms post-stimulus onset']);

if (~exist('include_reps', 'var')==1 | include_reps == 0)
	include_reps = 1:max_reps;
    if (use_testinfo & ~isempty(curtestinfo.include_reps))
	    include_reps = curtestinfo.include_reps;
    end;
end;
nincl_reps = length(include_reps);

% FILTER OUT UNWANTED SPIKES
if (isempty(data_array))
   disp(['Test ' num2str(test) ': Error processing data file.  Data array empty.']);
   mean_resp_surf = -1;
   return;
end;


% DO THE LATENCY SHIFT
if (~rc_flag & test_type~=300)
  data_array(:,1)-intrinsic_latency;
  data_array(:,1) = data_array(:,1)-intrinsic_latency;
  data_array = data_array(data_array(:,1)>=0,:);   % trim out negative values
end;

% CUT TO SPACE, IF APPROPRIATE

if (test_type==300)
   data_array(:,1) = data_array(:,1)-int_latency_space;
   data_array = data_array(data_array(:,1)>=0,:);   % trim out negative values
   [mean_resp_surf, std_surf, testpars] = ...
   		proc_space_test(test_head, tag, cluster, data_array, locs, start_count, end_count,...
   						max_reps, stim_dur, isi, rast_overrun, plot_flag, include_reps);
   dim2vals = locs;
   dim1vals = [];
   return; 

end;

% DRAW RASTER PLOT AND CREATE TITLE 
if (plot_flag==1 | plot_flag==2)
   figure
	rast_stim_dur = stim_dur;
	rast_total_width = stim_dur+rast_overrun;
	if (bg_flag)
		rast_stim_dur = rast_stim_dur + 2*abs(bg_lead);
		rast_total_width = rast_total_width + 2*abs(bg_lead);
	end;
	
	if (strcmp(dim1label, 'FREQ'))
		% cluge the data to make the freq axis plot from low to high
		temp = data_array(:,4);
		data_array(:,4) = max(data_array(:,4)) - data_array(:,4) + 1;
		ah = draw_grid(data_array(:,[1 3 4 5]), flipud(dim1vals), dim2vals, max_reps, rast_stim_dur, ...
			                 rast_total_width, dim1label, dim2label, count_starttime, count_endtime);	
		data_array(:,4) = temp;
	else
		ah = draw_grid(data_array(:,[1 3 4 5]), dim1vals, dim2vals, max_reps, rast_stim_dur, ...
			                 rast_total_width, dim1label, dim2label, count_starttime, count_endtime);	
	end;	
	
	fh = get(ah, 'parent');
	  
	n_rows = size(dim1vals,1);
	n_cols = size(dim2vals,1);
	rast_width = (isi + stim_dur)*1.05;	% see col_space in draw_grid.m
	rast_height = (max_reps+5) * 4;	 % numbers come from row_space and grid aspect ratio
												 % set in draw_grid.m
	if (n_cols == 1 | n_rows*rast_height > n_cols*rast_width)
		set(ah, 'position', [0.1 0 .9 .9]);
		set(fh, 'PaperOrientation', 'portrait');
		set(fh, 'PaperPosition', [.25 .25 8 10.5])
	else
		set(ah, 'position', [0.05 0 .95 .9]);
		% if title goes off page, try this: set(ah, 'Position', [0.07 0 1 .8])
		set(fh, 'PaperOrientation', 'landscape');
		set(fh, 'PaperPosition', [.25 .25 10.5 8])
	end;
	
	tempstr = [];
   tempstr = [tag '  SI ' num2str((2)) ...
				  '    T ' num2str(test) ...
				  '    Cl ' num2str(cluster)];
	if (test_head(1,4)==0) tempstr = [tempstr '   Tones '];end
   if(test_head(1,4)==1) tempstr=[tempstr '   Noise'];end
   if(test_head(1,4)==2) tempstr=[tempstr '   BP Noise'];end
	
	tempstr = [tempstr '    ISI=' num2str(isi) ...
				  '     STIMDUR=' num2str(stim_dur)];
	if (bg_flag)
		tempstr = [tempstr '  BG lead=' num2str(bg_lead)];
	end;
	if (spf_flag)
		tempstr = [tempstr '  Prefilt El ' num2str(spf_el) ', Az ' num2str(spf_az)]; 
	end;
	
	titlestr(1) = {tempstr};
	tempstr = [];
	if (size(ild_vals,1)==1) tempstr = ['ILD= ' num2str(ild_vals)]; end;
	if (size(itd_vals,1)==1) tempstr = [tempstr '     ITD= ' num2str(itd_vals)]; end;
	%if (size(phase_vals,1)==1) tempstr = [tempstr '	  PHASE= ' num2str(phase_vals)]; end;
	if (size(freq_vals,1)==1 & freq_vals(1)~=0) tempstr = [tempstr '	FREQ= ' num2str(freq_vals)]; end;
	if (size(abi_vals,1)==1) tempstr = [tempstr '     ABI= ' num2str(abi_vals)]; end;
	tempstr = [tempstr '       (ML: ' num2str(test_head(1, 6)) ',  AP: ' num2str(test_head(1,5)) ...
								',  Dpth: ' num2str(test_head(1,7)) ')'];
	titlestr(2) = {tempstr};
	txthndl = title(titlestr);
end; % if plot_flag

% COUNT SPIKES

disp('counting spikes...');

% DISPLAY NUMBER OF SPIKES FOR EACH REP 

for rc=1:max_reps
	rep_sum = sum(data_array(:,3)==rc);
	disp(['REP ' num2str(rc) ': ' num2str(rep_sum) ' spikes']);
end;

[mean_resp_surf std_surf] = count_spikes(data_array, dim1vals, dim2vals, start_count, end_count, include_reps);

% DRAW INSET RESP SURF PLOT

if (plot_flag==2 | plot_flag==3) 

  fh = gcf;
  if (plot_flag==2) 
     ah = axes;
  else
	 ah = gca
  end;
  if (size(dim2vals,1) == 1)
	  plot(dim1vals, mean_resp_surf);
	  ah = gca;
	  if (plot_flag==2) 
		  set(ah, 'position', [.65 .05 .3 .24]);
	  end;
	  set(ah, 'fontsize', 8);
	  cur_axis = axis;
	  axis([min(dim1vals) max(dim1vals) cur_axis(3) cur_axis(4)]);
  else
	  if (plot_flag ==2)
	  	set(ah, 'position', [.82 .02 .14 .18]);
	  end;
	  plotsurf(dim2vals, dim1vals, mean_resp_surf);
	  colorbar;
	  ah = gca;
	  set(ah, 'fontsize', 6);	 
	  if (strcmp(dim2label, 'ABI')) % if second axis (x axis) is ABI
		  set(ah, 'xdir', 'reverse');
	  end;
	  if (strcmp(dim1label, 'LEFT ATTEN'))
			%set(ah, 'ydir', 'normal');
			set(ah, 'xdir', 'reverse');
	  end;
  end;

  % CREATE THE RESP SURF TITLE

  titlestr = [];
  titlestr{1} = ['Counted: ' num2str(start_count) ...
			  ' to ' num2str(end_count) '	 '];
  num_inc_reps = size(include_reps, 2);
  if (num_inc_reps<=10) 
	  titlestr{1} = [titlestr{1} 'Reps: [' num2str(include_reps) ']']; 
  else
	  titlestr{1} = [titlestr{1} ...
			 'Reps: [' num2str(include_reps(1)) ' to ' num2str(include_reps(num_inc_reps)) ']']; 
  end;
  titlestr{2} = ['Spont Rate: ' num2str(spont_rate) ' spikes/stim'];
  if (plot_flag==2)
  	txthndl = title(titlestr);
  	set(txthndl, 'fontsize', 6);
  end;
  
 
  if (plot_flag==3)	% add title for stand-alone test
  
	  tempstr = [];
	  tempstr = [tag ' SI ' num2str(test_head(2)) ...
         ' T ' num2str(test) ...
         ' CL ' num2str(cluster)];
	  tempstr = [tempstr ' ISI=' num2str(isi) ...
				  ' DUR=' num2str(stim_dur)];
      if (spf_flag)
         tempstr = [tempstr '  Prefilt El ' num2str(spf_el) ', Az ' num2str(spf_az)]; 
      end;
	  
      titlestr2(1) = {tempstr};
	  tempstr = [];
	  if (size(ild_vals,1)==1) tempstr = ['ILD= ' num2str(ild_vals)]; end;
	  if (size(itd_vals,1)==1) tempstr = [tempstr '   ITD= ' num2str(itd_vals)]; end;
	  if (size(freq_vals,1)==1 & freq_vals(1)~=0) tempstr = [tempstr '   FREQ= ' num2str(freq_vals)]; end;
	  if (size(abi_vals,1)==1) tempstr = [tempstr '   ABI= ' num2str(abi_vals)]; end;
      if (test_head(1,4)==1) tempstr = [tempstr '   White Noise ']; end;
	  titlestr2(2) = {tempstr};

	  txthndl = title(titlestr2);
	  set(txthndl, 'fontsize', 8);
  end;

  xlabel(dim1label);
  if (size(dim2vals,1)==1)
	  xlabel(dim1label, 'fontsize', 8); 
	  ylabel('Mean Spikes/Trial','fontsize', 8); 
  else
	  xlabel(dim2label, 'fontsize', 6);
	  ylabel(dim1label, 'fontsize', 6);
  end;
  if (plot_flag==3);
     set(gca, 'fontsize', 7);
  end;
  
  %colormap(1-gray);
end;	% if plot_flag

% MAKE PARAMETER RETURN VALUES

if (size(freq_vals,1)==1) testpars(1)=freq_vals;
else testpars(1)=NaN;
end;

if (size(itd_vals, 1)==1) testpars(2)=itd_vals;
else testpars(2) = NaN;
end;

if (size(ild_vals, 1)==1) testpars(3)=ild_vals;
else testpars(3) = NaN;
end;

if (size(abi_vals, 1)==1) testpars(4)=abi_vals;
else testpars(4) = NaN;
end;

testpars(5:7)=test_head(5:7);

testpars(8) = stim_dur;



return;



function [mean_resp_surf, std_surf] = count_spikes(data_array, dim1vals, dim2vals, startcount, endcount, include_reps);
% COUNT_SPIKES
%	 takes a data array and counts spikes for each condition
%	 data array = [spike_time clust rep dim1val dim2val];

data_array = data_array(data_array(:,1)>startcount & data_array(:,1)<endcount, [3 4 5]);

for i=1:size(dim1vals,1);
	for j=1:size(dim2vals, 1);
		data_vals = data_array((data_array(:,2) == i & data_array(:,3) == j), 1);
		for rep_i = 1:size(include_reps,2);
			cur_rep = include_reps(rep_i);	
			reptab(rep_i) = sum(data_vals==cur_rep);
		end;
		number = size(reptab,1);
		mean_resp_surf(i, j) = mean(reptab);
		%stderr_surf(i, j) = std(reptab)/sqrt(number);
		std_surf(i, j) = std(reptab);    
	end;
end;


% ***************************
% PROC_SPACE_TEST
% ***************************

function [mean_resp_surf, std_surf, testpars] = ...
   		proc_space_test(header, tag, cluster, data_array, locs, start_count, end_count,...
   						max_reps, stim_dur, isi, rast_overrun, plot_flag, include_reps); 

testpars = [NaN NaN NaN header(4)];
testpars(5:7)=header(5:7);
testpars(8) = stim_dur;

% COUNT SPIKES

disp('counting spikes...');

if (~exist('include_reps', 'var')==1 | include_reps == 0)
	include_reps = 1:max_reps;
end;

% DISPLAY NUMBER OF SPIKES FOR EACH REP 

for rc=1:max_reps
	rep_sum = sum(data_array(:,3)==rc);
	disp(['REP ' num2str(rc) ': ' num2str(rep_sum) ' spikes']);
end;

% prepare data for raster
data_array_el0 = data_array(data_array(:,4)==0,:);
az_vals = data_array_el0(:,5);

az_axis = sort(az_vals);    % az_axis is single values of all azimuth positions
az_i = abs([1; diff(az_axis)])>0;
az_axis = az_axis(az_i);

az_vals_i = zeros(size(az_vals));
for i = 1:size(az_axis)
   az_vals_i(az_vals==az_axis(i)) = i*ones(size(az_vals(az_vals==az_axis(i))));
end;

% convert data_array_el0 to standard format where 4th and 5th cols carry
% ordinal condition data for x and y axes
data_array_el0(:,4) = az_vals_i;  % converts from raw values to index into those values
data_array_el0(:,5) = ones(size(data_array_el0,1),1); 

% DRAW THE RASTER FOR AZIMUTH ZERO

if (plot_flag==1 | plot_flag==2)
  figure
  rast_stim_dur = stim_dur;
  rast_total_width = stim_dur+rast_overrun;

  ah = draw_grid(data_array_el0(:,[1 3 4 5]), az_axis, 0, max_reps, rast_stim_dur, ...
                 rast_total_width, 'azimuth', '', start_count, end_count);	
  fh = get(ah, 'parent');
  
  %n_rows = size(dim1vals,1);
  %n_cols = 1;
  rast_width = (isi + stim_dur)*1.05;	% see col_space in draw_grid.m
  rast_height = (max_reps+5) * 4;	 % numbers come from row_space and grid aspect ratio
											 % set in draw_grid.m
  set(ah, 'position', [0.1 0 .9 .9]);
  set(fh, 'PaperOrientation', 'portrait');
  set(fh, 'PaperPosition', [.25 .25 8 10.5]);

  tempstr = [];
  tempstr = [tag '    SI ' num2str(header(2)) ...
			  '    T ' num2str(header(1)) ...
			  '    Cl ' num2str(cluster)];
  tempstr = [tempstr '   ISI=' num2str(isi) ...
			  '    STIMDUR=' num2str(stim_dur)];
  % add info about test type (true, flat itd or abi)
  if (size(header,2)==13)
     if (header(11)==0 & header(12)==0) tempstr = [tempstr '  true space']; end;
     if (header(11)==1 & header(12)==0) tempstr = [tempstr '  ABI equal']; end;
     if (header(11)==0 & header(12)==1) tempstr = [tempstr '  ITD equal']; end;
     if (header(11)==1 & header(12)==1) tempstr = [tempstr '  ABI & ITD equal']; end;
  end;			  
  titlestr(1) = {tempstr};
  tempstr = [];
  tempstr = [tempstr '  ABI= ' num2str(header(4))];
  if (size(header,2)==13 & header(12)==1)
     tempstr = [tempstr '   ITD= ' num2str(header(13))];
  end;
  tempstr = [tempstr '     (ML: ' num2str(header(1, 6)) ',  AP: ' num2str(header(1,5)) ...
							',  Dpth: ' num2str(header(1,7)) ')'];

  titlestr(2) = {tempstr};
  txthndl = title(titlestr);
end; % if plot_flag says "rasters"

data_array = data_array(data_array(:,1)>start_count & data_array(:,1)<end_count, [3 4 5]);

for i=1:size(locs,1);
	data_vals = data_array((data_array(:,2) == locs(i,1) & data_array(:,3) == locs(i,2)), 1);
	for rep_i = 1:size(include_reps,2);
		cur_rep = include_reps(rep_i);
		reptab(rep_i) = sum(data_vals==cur_rep);
	end;
	number = size(reptab,1);
	mean_resp_surf(i) = mean(reptab);
	%stderr_surf(i) = std(reptab)/sqrt(number);
    std_surf(i) = std(reptab);
   
end;

if (plot_flag==2 | plot_flag==3) 

	fh = gcf;
	if (plot_flag==2) 
	   ah = axes;
	else
	   ah = gca;
	end;
   
	map = locs';
	
	[azi, ele, data] = array2diamond(mean_resp_surf, map);
	
	[AZ EL] = meshgrid(azi, ele);
	
	% interpolate missing values
	
	% generate mask for missing points
	missmask = NaN*ones(size(data));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(locs(:,1)==el & locs(:,2)==az))
	       missmask(AZ==az & EL==el) = 1;
	    else
	       missmask(AZ==az & EL==el) = 0;
	    end; 
	  end;
	end;
	
	% replace all valid locations presently containing NaN's with zeros
	ri = (missmask==1);
	data(ri) = zeros(size(data(ri)));
	
	% generate interpolation function
	intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
	intval = conv2(data, intfun,'same')/4;
	intval = intval.*missmask;
	intval(find(isnan(intval))) = zeros(size(find(isnan(intval))));
	
	newdata = intval + data;
	
	% display parts not shown in using pcolor and adjust axes
	elstep = min(diff(ele));
	azstep = min(diff(azi));
	n_azi = size(newdata,2);
	n_ele = size(newdata,1);
	disp_el = [ele ele(n_ele) + elstep];
	disp_az = [azi  azi(n_azi)  + azstep];
	disp_data = [newdata newdata(:,n_azi)];
	disp_data = [disp_data; disp_data(n_ele,:)];
	
	pcolor(disp_az,disp_el,disp_data);
	shading flat;
	cbh = colorbar;
	whitebg([0 0 0]);
	
	% realign axes
	yaxis_shift = elstep/2;
	xaxis_shift = azstep/2;
	ca = gca;
	xtl = get(ca,'xticklabel');
	ytl = get(ca,'yticklabel');
	set(ca,'xtick',get(ca,'xtick')+xaxis_shift);
	set(ca,'ytick',get(ca,'ytick')+yaxis_shift);
	set(ca,'xticklabel',xtl);
	set(ca,'yticklabel',ytl);
	set(ca,'dataaspectratio', [1 1 1])
	
	% position and font size for graph and colorbar
	if (plot_flag==2) 
	   set(ah, 'fontsize', 8);
	   set(ah, 'position', [.65 .05 .3 .3]);
	   set(cbh, 'fontsize', 8);
	   set(cbh, 'position', [.96 .096 .02 .207]);
	end;
	
	
	% make the title
	
	titlestr = [];
	titlestr = ['Counted: ' num2str(start_count) ...
				  ' to ' num2str(end_count) '  '];
	num_inc_reps = size(include_reps, 2);
	if (num_inc_reps<=10) 
		  titlestr = [titlestr 'Reps: [' num2str(include_reps) ']']; 
	else
		  titlestr = [titlestr ...
				 'Reps: [' num2str(include_reps(1)) ' to ' num2str(include_reps(num_inc_reps)) ']']; 
	end;

    % add more info for stand-alone test
	if (plot_flag~=3)
	
		txthndl = title(titlestr);
		set(txthndl, 'fontsize', 6);
	    tempstr = [tag '  Site ' num2str(header(2)) ...
					  '   Test ' num2str(header(1)) ...
					  '   Clust ' num2str(cluster)];
		tempstr=[tempstr '   Noise'];
		tempstr = [tempstr '  ISI=' num2str(isi) ...
			  '    STIMDUR=' num2str(stim_dur)];
		titlestr2(1) = {tempstr};
		titlestr2(2) = {titlestr};
		tempstr = [];
		txthndl = title(titlestr2);
	else
		tempstr = [];
		tempstr = [tag ' SI ' num2str(header(2)) ...
	         ' T ' num2str(header(1)) ...
	         ' CL ' num2str(cluster)];
		tempstr = [tempstr ' ISI=' num2str(isi) ...
					  ' DUR=' num2str(stim_dur)];
		  
	    titlestr2(1) = {tempstr};
		tempstr = [];
		tempstr = [tempstr 'ABI= ' num2str(header(4))];
        if (size(header,2)==13 & header(12)==1)
           tempstr = [tempstr '   ITD= ' num2str(header(13))];
        end;
        if (size(header,2)==13)
	        if (header(11)==0 & header(12)==0) tempstr = [tempstr '  true space']; end;
	        if (header(11)==1 & header(12)==0) tempstr = [tempstr '  ABI equal']; end;
	        if (header(11)==0 & header(12)==1) tempstr = [tempstr '  ITD equal']; end;
	        if (header(11)==1 & header(12)==1) tempstr = [tempstr '  ABI & ITD equal']; end;
        end;		
		titlestr2(2) = {tempstr};
	
		txthndl = title(titlestr2);
		set(txthndl, 'fontsize', 6);
	end;
	
	
	if (plot_flag==3);
	    set(gca, 'fontsize', 7);
		xlabel('azimuth', 'fontsize', 7); 
		ylabel('elevation', 'fontsize', 7);
		
	else
		xlabel('azimuth'); 
		ylabel('elevation');
	end;
	
	if (plot_flag==3) set(gcf, 'paperposition', [2.3 3.47 3.9 4]); end;
	%colormap(1-gray);

end; % if inset plot


% ************************** END ******************************


function y = proc_raw(rec_filename, test);
% PROC_RAW	 process raw test data and store the results in a cache file
%[resp_surf, stderr_surf, xaxis, yaxis] = proc_raw(rec_filename, test);
%	This version stashes the loaded data for each analyzed test in a file so that
%	the next time that test is analyzed the load will be very fast.

% note: we can use set(1, 'paperposition', [3 .5 3 10]) to set the raster display
% to display better [left bottom width height]

% set data directory

datadirectory = 'd:\mlspezio\physioldata\899\';
cutdirectory = '899cut';
cache_files_directory = 'd:\mlspezio\matlab\scripts\physiol_analysis\899analysis\899_cache';

disp(['loading rec file '[datadirectory rec_filename] '...']);
[b, header, files, cond_data, trial_order, spikes, comments] = readrec([datadirectory rec_filename], test);

% check that required condition data is in this rec file
if (b.number==-1)
   disp(['Test ' int2str(test) ' not found in specified record file.']);
   y = -1;
   return;
end;
disp('done.');

% set some constants

test_type = header(3);

if (test_type>1100)
  disp(['Unable to process test ' int2str(test) '...aborting.']);
  y = -1;
  return;
end;

test_class = fix(test_type/100);
test_type = test_type - test_class*100;

num_trials = size(cond_data,1);
isi = header(9);
stim_dur = header(10);
max_reps = max(cond_data(:,2));
tag = sprintf('%.0f%c', b.number, b.side);
site = header(2);

if (test_class==3)
   disp('processing space test...');
   y = proc_space_raw(b, header, files, cond_data, spikes, comments, ...
                      datadirectory, cutdirectory, cache_files_directory); 
   return;
end;

rc_flag = 0;
if (test_class==1) 
	rc_flag = 1;
	rast_overrun = 20;
	num_trials = size(trial_order,2);
	trials_per_block = size(cond_data,1);
	max_reps = num_trials/trials_per_block;
end;

bg_flag = 0;
bg_lead = 0;
if (test_class == 2)
   bg_flag = 1;
   bg_az = header(12);
   bg_el = header(13);
   bg_abi = header(14);
   bg_lead = header(15);
end;

spf_flag = 0;
if (test_class == 10)
   spf_flag = 1;
   spf_el = header(12);
   spf_az = header(13);
end;

% find bird letter for cut file name
if (b.number==874) bird_letter='C'; end;
if (b.number==885) bird_letter='D'; end;
if (b.number==855) bird_letter='E'; end;
if (b.number==890) bird_letter='F'; end;
if (b.number==891) bird_letter='G'; end;
if (b.number==897) bird_letter='H'; end;

% load triggers and test header/trial condition data

test_head = header;
rep1_cond = cond_data(cond_data(:,2)==1, :);

if (size(cond_data, 1)<2)
	disp('Not enough condition data to process test.  Aborting...');
	y = -1;
	return;
end;

dw_ticks_per_ms = 10;
dw_tdt_cor = 1.00158;

% Select all levels of every variable, sort and eliminate redundancy

ild_vals = sort(rep1_cond(:,3));
ild_i =abs([1; diff(ild_vals)])>0;
ild_vals = ild_vals(ild_i);

itd_vals = sort(rep1_cond(:,4));
itd_i =abs([1; diff(itd_vals)])>0;
itd_vals = itd_vals(itd_i);

phase_vals = sort(rep1_cond(:,5));
phase_i =abs([1; diff(phase_vals)])>0;
phase_vals = phase_vals(phase_i);

freq_vals = sort(rep1_cond(:,6));
freq_i =abs([1; diff(freq_vals)])>0;
freq_vals = freq_vals(freq_i);

abi_vals = sort(rep1_cond(:,7));
abi_i =abs([1; diff(abi_vals)])>0;
abi_vals = abi_vals(abi_i);
abi_vals = flipud(abi_vals);

latt_vals = sort(rep1_cond(:,8));
latt_i =abs([1; diff(latt_vals)])>0;
latt_vals = latt_vals(latt_i);

ratt_vals = sort(rep1_cond(:,9));
ratt_i =abs([1; diff(ratt_vals)])>0;
ratt_vals = ratt_vals(ratt_i);
ratt_vals = flipud(ratt_vals); % insures rasters plotted in order of increasing loudness

% test types: 1=ILD, 2=ITD, 3=PHASE, 4=FREQ, 5=ABI, 6=ILD-ITD
% 7 = ILD-FREQ, 8 = ILD-ABI, 9 = ITD-FREQ, 10 = FREQ-ABI, 11 = PHASE-ABI, 12 = RATT-LATT,
% 13 = ILD-ITD-FREQ


% FORM THE DATA ARRAY
% Match all spikes to their corresponding condition values
% and subtract off time at beginning of trial so that times start from zero
% data_array = [spike_time rep dim1valnum dim2valnum]

if (test_type==1)	 % ILD
  cond_data_col1 = 3;
  cond_data_col2 = 0;
  dim1vals = ild_vals;
  dim2vals = 0;
  dim1label = 'ILD';
  dim2label = 'null';
end;

if (test_type==2)	 % ITD
  cond_data_col1 = 4;
  cond_data_col2 = 0;
  dim1vals = itd_vals;
  dim2vals = 0;
  dim1label = 'ITD';
  dim2label = 'null';
end;

if (test_type==4)	 % FREQ
  cond_data_col1 = 6;
  cond_data_col2 = 0;
  dim1vals = freq_vals;
  dim2vals = 0;
  dim1label = 'FREQ';
  dim2label = 'null';
end;

if (test_type==5)	 % ABI
  cond_data_col1 = 7;
  cond_data_col2 = 0;
  dim1vals = abi_vals;
  dim2vals = 0;
  dim1label = 'ABI'; 
  dim2label = 'null';
end;

if (test_type==6)	 % ILD-ITD
  cond_data_col1 = 4;
  cond_data_col2 = 3;
  dim1vals = itd_vals;
  dim2vals = ild_vals;
  dim1label = 'ITD';
  dim2label = 'ILD';
end;

if (test_type==7)	 % ILD-FREQ
  cond_data_col1 = 6;
  cond_data_col2 = 3;
  dim1vals = freq_vals;
  dim2vals = ild_vals;
  dim1label = 'FREQ';
  dim2label = 'ILD';
end;

if (test_type==9) % ITD-FREQ
  cond_data_col1 = 4;
  cond_data_col2 = 6;
  dim1vals = itd_vals;
  dim2vals = freq_vals;
  dim1label = 'ITD';
  dim2label = 'FREQ';
end;

if (test_type==10) % FREQ-ABI
  cond_data_col1 = 6;
  cond_data_col2 = 7;
  dim1vals = freq_vals;
  dim2vals = abi_vals;
  dim1label = 'FREQ';
  dim2label = 'ABI';	 % this value is used later to see whether to invert the x axis
end;

if (test_type==12) % LATT-RATT
  cond_data_col1 = 8;
  cond_data_col2 = 9;
  dim1vals = latt_vals;
  dim2vals = ratt_vals;
  dim1label = 'LEFT ATTEN';
  dim2label = 'RIGHT ATTEN';
end;

% attempt to load cut data and if not load from spike file

disp('loading data file...');
[trigger clust0 clust1 clust2 clust3 clust4] = ...
readact(test, bird_letter, b.side, site, cutdirectory); 

cut_data_flag = 1;
if (trigger == -1) 
	  
   cut_data_flag = 0;
   disp('cut file not available...using spike file');
   % get spike data from spike file
   % first need to strip off original pathname
   spikefile = files{3};
   for i=size(spikefile,2):-1:1; 
      if (spikefile(i)=='\') break; end;
   end;
   
   spikefile = [datadirectory '\' spikefile(i+1:size(spikefile,2))]

   clust1 = spikein(spikefile)';
   trigger = 0;
   if (clust1==-1)
      disp(['Error reading spike data file for test ' int2str(test) '.  Aborting...']);
      y = -1;
      return;
   end;
   % convert for MII/TDT timing mismatch
   clust1 = clust1*1.00166;
   clust1 = clust1/1000;  % do it in ms
  
	clust0 = [];
	clust2 = [];
	clust3 = [];
	clust4 = [];

end;

disp('done.');

clust_count(1) = size(clust0, 2);
clust_count(2) = size(clust1, 2);
clust_count(3) = size(clust2, 2);
clust_count(4) = size(clust3, 2);
clust_count(5) = size(clust4, 2);

disp(['Tabulating spikes...']);
data_array = [];

% PROCESS CUT FILE STANDARD AND BACKGROUND NOISE
if (~rc_flag & cut_data_flag)

	trigger = trigger';
	
	if (size(trigger,1) ~= size(cond_data,1))
		disp('***** Warning:	Trial conditions and triggers dont match!');
		disp('trial data:');
		size(cond_data,1)
		disp('triggers:');
		size(trigger,1)
		datasize = min(size(trigger,1), size(cond_data,1));
		trigger = trigger(1:datasize, :);
		cond_data = cond_data(1:datasize, :);
		num_trials = size(cond_data,1);
	end;

	for tr=1:num_trials;
   	d1n = find(dim1vals == cond_data(tr,cond_data_col1));
		if (dim2vals==0) 
			d2n = 1;
		else 
	  		d2n = find(dim2vals == cond_data(tr,cond_data_col2));
		end;
		rep = cond_data(tr,2);
		trigger1 = trigger(tr,1);
		if (tr+1<= num_trials) 
			trigger2 = trigger(tr+1,1);
		else
			trigger2 = trigger1+stim_dur+isi+2*abs(bg_lead)+10;
		end;
		for cl = [0:4]
			spike_var = ['clust' int2str(cl)];
			eval(['spikes = ' spike_var ';']);
			spikes = spikes';
			if (~isempty(spikes))
				trial_spikes = spikes(spikes>trigger1 & spikes<trigger2);
				trial_spikes = (trial_spikes - trigger1)/dw_ticks_per_ms;
				if (trial_spikes)
					data_array = [data_array; ...
						trial_spikes ones(size(trial_spikes))*[cl rep d1n d2n]];
				end; % if trial_spikes
			end; % if ~empty(spikes)
		end; % for cl			
  	end;  % for tr
  	
end;

% PROCESS SPIKE FILE STANDARD OR BACKGROUND NOISE
if (~rc_flag & ~cut_data_flag)
	
	trial_starts = find(clust1==0)';
	cl = 1;
	
	if (size(trial_starts,1) ~= size(cond_data,1))
		disp('***** Warning:	Trial conditions and triggers dont match!');
		disp('trial data:');
		size(cond_data,1)
		disp('triggers:');
		size(trial_starts,1)
		datasize = min(size(trial_starts,1), size(cond_data,1));
		trial_starts = trial_starts(1:datasize);
		cond_data = cond_data(1:datasize, :);
		num_trials = size(cond_data,1);
	end;
	for tr=1:num_trials;
   	d1n = find(dim1vals == cond_data(tr,cond_data_col1));
		if (dim2vals==0) 
			d2n = 1;
		else 
	  		d2n = find(dim2vals == cond_data(tr,cond_data_col2));
		end;
		rep = cond_data(tr,2);
		if (tr==num_trials)
			trial_end = size(trial_starts,1);
		else
			trial_end = trial_starts(tr+1)-1;
		end;
      trial_spikes = clust1(trial_starts(tr)+1:trial_end)';
      
	   if (trial_spikes)
			data_array = [data_array; ...
						trial_spikes ones(size(trial_spikes))*[cl rep d1n d2n]];
		end; % if trial_spikes
   end;  % for tr
  	
end;  % process spike file


if (cut_data_flag)  % convert to ms
   trigger = dw_tdt_cor*trigger/dw_ticks_per_ms;
   clust0 = dw_tdt_cor*clust0/dw_ticks_per_ms;
   clust1 = dw_tdt_cor*clust1/dw_ticks_per_ms;
	clust2 = dw_tdt_cor*clust2/dw_ticks_per_ms;
	clust3 = dw_tdt_cor*clust3/dw_ticks_per_ms;
	clust4 = dw_tdt_cor*clust4/dw_ticks_per_ms;
end;

% PROCESS REV COR FROM CUT OR SPIKE FILE
if (rc_flag)
   
   % check for abort
  
   for i=1:size(comments,2)
      if (isempty(comments{1})) break; end;
      if (strcmp(comments{i}(1:3),'Rep'))
         disp(comments);
         actual_reps = sscanf(comments{i}(7:size(comments{i},2)),'%d',1);
         num_trials = trials_per_block*(actual_reps-1);
         max_reps = actual_reps-1;
         if (num_trials<=0)
           disp('test was aborted...terminating processing.');
           return;
         end;
   	end;
   end;

	onset = (1:size(trial_order,2))*stim_dur;
	cl = 0;
	%lag = stim_dur-22;
	lag = 89
	test_onset = min(trigger);
	clust0 = clust0 - test_onset + lag;
	clust1 = clust1 - test_onset + lag;
	clust2 = clust2 - test_onset + lag;
	clust3 = clust3 - test_onset + lag;
	clust4 = clust4 - test_onset + lag;
	
	for tr=1:num_trials;
		d1n = find(dim1vals == cond_data(trial_order(tr),cond_data_col1));
		if (dim2vals==0) 
		  d2n = 1;
		else 
		  d2n = find(dim2vals == cond_data(trial_order(tr),cond_data_col2));
		end;
		rep = floor((tr-1)/trials_per_block)+1;
		offset = onset(tr)+stim_dur+rast_overrun;
	   for cl = [0:4]
			spike_var = ['clust' int2str(cl)];
			eval(['spikes = ' spike_var ';']);
			spikes = spikes';
			if (~isempty(spikes)) 
	      	trial_spikes = spikes(spikes>onset(tr) & spikes<offset);
	      	trial_spikes = (trial_spikes - onset(tr));
	      	if (trial_spikes)
	         	data_array = [data_array; ...
	         					trial_spikes ones(size(trial_spikes))*[cl rep d1n d2n]];
				end; % if trial_spike_data
			end; % if spikes
		end;  % for cl
	end; % for tr
		
end;

% round frequencies	(note this must be done after processing spikes because
%                     processing relies on matching dim vals to entries in record		 
if (cond_data_col1==6) 
	dim1vals = round(dim1vals); 
end;
if (cond_data_col2==6) 
	dim2vals = round(dim2vals);
end;
						
% TABULATE SPONTANEOUS SPIKES BEFORE AND AFTER TEST

beforetime_i = min([min(clust0) min(clust1) min(clust2) min(clust3) min(clust4)]);
beforetime_f = min(trigger);
if (rc_flag)
	aftertime_i = onset(size(onset,2)) + stim_dur;
else
	aftertime_i = max(trigger) + (stim_dur+isi+2*abs(bg_lead));
end;
aftertime_f = max([max(clust0) max(clust1) max(clust2) max(clust3) max(clust4)]);
beforedur = (beforetime_f - beforetime_i);
afterdur = (aftertime_f - aftertime_i);
if (beforedur<0) beforedur = 0; end;
if (afterdur<0) afterdur = 0; end;
		
prespikes = [];
postspikes = [];
for cl = [0:4]
	spike_var = ['clust' int2str(cl)];
	eval(['spikes = ' spike_var ';']);
	spikes = spikes';
	if (~isempty(spikes))
		
		beforespikes = spikes(spikes>beforetime_i & spikes<beforetime_f);
		beforespikes = beforespikes - beforetime_i;
		afterspikes = spikes(spikes>aftertime_i & spikes<aftertime_f);
		afterspikes = afterspikes - aftertime_i;
		
		prespikes = [prespikes; beforespikes ones(size(beforespikes))*cl];
		postspikes = [postspikes; afterspikes ones(size(afterspikes))*cl];
		
	end; % if ~empty(spikes)
end; % for cl

% TRIM OUT THE ZERO CLUSTER BECAUSE IT TAKES UP SPACE AND I NEVER USE IT
keep_index = ~data_array(:,2)==0;
data_array = data_array(keep_index,:);


% SAVE THE RESULTS TO CACHE FILE

data_file = sprintf([tag 'SITE%02gT%03g.MAT'], site, test);
disp(['Saving file: ' data_file]);
data_file = [cache_files_directory '\' data_file];
test_head = header;

eval(['save ' data_file ' test_head test_type data_array ' ...
                        ' dim2vals dim1vals files' ...
					    ' max_reps dim1label dim2label' ...
						' tag stim_dur isi clust_count' ...
						' prespikes postspikes beforedur afterdur '...
						' rc_flag bg_flag spf_flag' ...
						' ild_vals itd_vals phase_vals freq_vals abi_vals;']);


if (bg_flag)
   eval(['save ' data_file ' bg_az bg_el bg_abi bg_lead -append;']);
end;

if (spf_flag)
   eval(['save ' data_file ' spf_az spf_el -append;']);
end;

return;

% ********************************************
% PROC_SPACE_RAW    
% ********************************************

function y = proc_space_raw(b, header, files, cond_data, spikes, comments,...
 									 datadirectory, cutdirectory, cache_files_directory); 
% proc_space_raw  processes a raw space file

rc_flag = 0;
test = header(1);
num_trials = size(cond_data,1);
isi = header(8);

stim_dur = header(9);
max_reps = max(cond_data(:,2));

tag = sprintf('%.0f%c', b.number, b.side);
site = header(2);
test_type = header(3);

% find bird letter for cut file name
if (b.number==874) bird_letter='C'; end;
if (b.number==885) bird_letter='D'; end;
if (b.number==855) bird_letter='E'; end;
if (b.number==890) bird_letter='F'; end;
if (b.number==891) bird_letter='G'; end;
if (b.number==897) bird_letter='H'; end;

% load triggers and test header/trial condition data

test_head = header;

if (size(cond_data, 1)<2)
	disp('Not enough condition data to process test.  Aborting...');
	y = -1;
	return;
end;

dw_ticks_per_ms = 10;
dw_tdt_cor = 1.00158;

% get the locations and sort
locs = cond_data(cond_data(:,2)==1,3:4);  % select first block

[sorted_array, loci] = sort(locs(:,2));
locs = locs(loci,:);
[sorted_array, loci] = sort(locs(:,1));
locs = locs(loci,:);

% load spikes from ACT or SPI file

disp('loading data file...');
[trigger clust0 clust1 clust2 clust3 clust4] = ...
readact(test, bird_letter, b.side, site, cutdirectory); 

cut_data_flag = 1;
if (trigger == -1) 
	  
   cut_data_flag = 0;
   disp('cut file not available...using spike file');
   % get spike data from spike file
   % first need to strip off original pathname
   spikefile = files{4};
   for i=size(spikefile,2):-1:1; 
      if (spikefile(i)=='\') break; end;
   end;
   
   spikefile = [datadirectory '\' spikefile(i+1:size(spikefile,2))];

   clust1 = spikein(spikefile)';
   trigger = 0;
   if (clust1==-1)
      disp(['Error reading spike data file for test ' int2str(test) '.  Aborting...']);
      y = -1;
      return;
   end;
   % convert for MII/TDT timing mismatch
   clust1 = clust1*1.00166;
   clust1 = clust1/1000;  % do it in ms
  
	clust0 = [];
	clust2 = [];
	clust3 = [];
	clust4 = [];

end;

clust_count(1) = size(clust0, 2);
clust_count(2) = size(clust1, 2);
clust_count(3) = size(clust2, 2);
clust_count(4) = size(clust3, 2);
clust_count(5) = size(clust4, 2);

disp(['Tabulating spikes...']);
data_array = [];

% PROCESS CUT FILE STANDARD AND BACKGROUND NOISE
if (cut_data_flag)

	trigger = trigger';
	
	if (size(trigger,1) ~= size(cond_data,1))
		disp('***** Warning:	Trial conditions and triggers dont match!');
		disp('trial data:');
		size(cond_data,1)
		disp('triggers:');
		size(trigger,1)
		datasize = min(size(trigger,1), size(cond_data,1));
		trigger = trigger(1:datasize, :);
		cond_data = cond_data(1:datasize, :);
		num_trials = size(cond_data,1);
	end;

	for tr=1:num_trials;
	   trialel = cond_data(tr,3);
	   trialaz = cond_data(tr,4);
		rep = cond_data(tr,2);
		trigger1 = trigger(tr,1);
		if (tr+1<= num_trials) 
			trigger2 = trigger(tr+1,1);
		else
			trigger2 = trigger1+stim_dur+isi+10;
		end;
		for cl = [0:4]
			spike_var = ['clust' int2str(cl)];
			eval(['spikes = ' spike_var ';']);
			spikes = spikes';
			if (~isempty(spikes))
				trial_spikes = spikes(spikes>trigger1 & spikes<trigger2);
				trial_spikes = (trial_spikes - trigger1)/dw_ticks_per_ms;
				if (trial_spikes)
					data_array = [data_array; ...
						trial_spikes ones(size(trial_spikes))*[cl rep trialel trialaz]];
				end; % if trial_spikes
			end; % if ~empty(spikes)
		end; % for cl			
  	end;  % for tr  	
  	
end;


% PROCESS SPIKE FILE STANDARD OR BACKGROUND NOISE
if (~cut_data_flag)
	
	trial_starts = find(clust1==0)';
	cl = 1;
	
	if (size(trial_starts,1) ~= size(cond_data,1))
		disp('***** Warning:	Trial conditions and triggers dont match!');
		disp('trial data:');
		size(cond_data,1)
		disp('triggers:');
		size(trial_starts,1)
		datasize = min(size(trial_starts,1), size(cond_data,1));
		trial_starts = trial_starts(1:datasize);
		cond_data = cond_data(1:datasize, :);
		num_trials = size(cond_data,1);
	end;
	for tr=1:num_trials;
   	
   	trialel = cond_data(tr,3);
	   trialaz = cond_data(tr,4);
		rep = cond_data(tr,2);
		
		if (tr==num_trials)
			trial_end = size(trial_starts,1);
		else
			trial_end = trial_starts(tr+1)-1;
		end;
      trial_spikes = clust1(trial_starts(tr)+1:trial_end)';
		
		if (trial_spikes)
					data_array = [data_array; ...
						trial_spikes ones(size(trial_spikes))*[cl rep trialel trialaz]];
		end; % if trial_spikes
		
   end;  % for tr
  	
end;  % process spike file

% TABULATE SPONTANEOUS SPIKES BEFORE AND AFTER TEST

beforetime_i = min([min(clust0) min(clust1) min(clust2) min(clust3) min(clust4)]);
beforetime_f = min(trigger);
if (rc_flag)
	aftertime_i = onset(size(onset,2)) + stim_dur;
else
	aftertime_i = max(trigger) + (stim_dur+isi);  % note computations done in dw ticks
end;
aftertime_f = max([max(clust0) max(clust1) max(clust2) max(clust3) max(clust4)]);
beforedur = (beforetime_f - beforetime_i)/dw_ticks_per_ms; % now in ms
afterdur = (aftertime_f - aftertime_i)/dw_ticks_per_ms;
if (beforedur<0) beforedur = 0; end;
if (afterdur<0) afterdur = 0; end;
		
prespikes = [];
postspikes = [];
for cl = [0:4]
	spike_var = ['clust' int2str(cl)];
	eval(['spikes = ' spike_var ';']);
	spikes = spikes';
	if (~isempty(spikes))
		
		beforespikes = spikes(spikes>beforetime_i & spikes<beforetime_f);
		beforespikes = beforespikes - beforetime_i;
		afterspikes = spikes(spikes>aftertime_i & spikes<aftertime_f);
		afterspikes = afterspikes - aftertime_i;
		
		prespikes = [prespikes; beforespikes ones(size(beforespikes))*cl];
		postspikes = [postspikes; afterspikes ones(size(afterspikes))*cl];
		
	end; % if ~empty(spikes)
end; % for cl


% TRIM OUT THE ZERO CLUSTER BECAUSE IT TAKES UP SPACE AND I NEVER USE IT
keep_index = ~data_array(:,2)==0;
data_array = data_array(keep_index,:);


% SAVE THE RESULTS TO CACHE FILE

data_file = sprintf([tag 'SITE%02gT%03g.MAT'], site, test);
disp(['Saving file: ' data_file]);
data_file = [cache_files_directory '\' data_file];
test_head = header;
eval(['save ' data_file ' test_head test_type data_array ' ...
                        ' locs files ' ...
					    ' max_reps ' ...
						' tag stim_dur isi clust_count' ...
						' prespikes postspikes beforedur afterdur '...
						' rc_flag; ']);

y = 1;  % return ok

return;
